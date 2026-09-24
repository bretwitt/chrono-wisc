// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2026 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: bgwitt
// =============================================================================

#include "chrono_planet/dust/ChDustField.h"

#include <algorithm>
#include <cmath>
#include <fstream>

#include "chrono/core/ChVector2.h"

#include "chrono_planet/core/Parallel.h"

namespace chrono {
namespace planet {

namespace {

// Height cache entries kept before the cache is cleared
constexpr size_t kMaxCachedHeights = 4000000;

// A sweep toward the Sun: slices across the dominant axis a of the direction, each voxel reading the point one
// slice upstream, which lies off its own row by (off_b, off_c) voxels along the other two axes.
struct Sweep {
    int a, b, c;
    int step;             // upstream slice is ia + step
    double off_b, off_c;  // offset of the upstream point in the other two axes (voxels)
    double length;        // distance to the upstream point (m)
};

Sweep MakeSweep(const ChVector3d& L, double voxel) {
    const double v[3] = {L.x(), L.y(), L.z()};
    int a = 0;
    for (int k = 1; k < 3; ++k)
        if (std::abs(v[k]) > std::abs(v[a]))
            a = k;
    Sweep s;
    s.a = a;
    s.b = (a + 1) % 3;
    s.c = (a + 2) % 3;
    s.step = v[a] > 0 ? 1 : -1;
    const double t = 1.0 / std::abs(v[a]);  // along L, in voxels, to move one voxel along a
    s.off_b = v[s.b] * t;
    s.off_c = v[s.c] * t;
    s.length = t * voxel;
    return s;
}

}  // namespace

ChDustField::ChDustField(const Params& params, HeightFunction ground) : m_params(params), m_ground(std::move(ground)), m_rng(params.seed) {
    double total = 0;
    for (size_t b = 0; b < m_params.bins.size(); ++b) {
        m_mass_extinction.push_back(3.0 * m_params.extinction_efficiency / (4.0 * m_params.grain_density * m_params.bins[b].radius));
        m_bin_share.push_back(m_params.bins[b].mass_fraction * m_mass_extinction[b]);
        total += m_bin_share[b];
    }
    for (double& share : m_bin_share)
        share = total > 0 ? share / total : 0.0;
}

double ChDustField::MassExtinction(int bin) const {
    return m_mass_extinction[bin];
}

ChVector3d ChDustField::Position(const Particle& p, double time) const {
    const double t = time - p.time0;
    return p.pos0 + p.vel0 * t + ChVector3d(0, 0, -0.5 * m_params.gravity * t * t);
}

double ChDustField::LatticeHeight(int i, int j) {
    const std::int64_t key = (static_cast<std::int64_t>(i) << 32) | static_cast<std::uint32_t>(j);
    auto it = m_heights.find(key);
    if (it != m_heights.end())
        return it->second;
    if (m_heights.size() >= kMaxCachedHeights)
        m_heights.clear();
    const float h = static_cast<float>(m_ground(i * m_lattice, j * m_lattice));
    m_heights.emplace(key, h);
    return h;
}

double ChDustField::Ground(double x, double y) {
    if (m_lattice != m_spec.voxel) {
        m_lattice = m_spec.voxel;
        m_heights.clear();
    }
    const double fx = x / m_lattice;
    const double fy = y / m_lattice;
    const int i = static_cast<int>(std::floor(fx));
    const int j = static_cast<int>(std::floor(fy));
    const double u = fx - i;
    const double v = fy - j;
    return (1 - u) * (1 - v) * LatticeHeight(i, j) + u * (1 - v) * LatticeHeight(i + 1, j) + (1 - u) * v * LatticeHeight(i, j + 1) + u * v * LatticeHeight(i + 1, j + 1);
}

void ChDustField::Emit(const ChVector3d& pos, const ChVector3d& vel, double mass, int bin, double time) {
    // The budget counts the mass as stored, so it balances exactly
    const float stored = static_cast<float>(mass);
    if (stored <= 0)
        return;
    m_stats.emitted += stored;
    if (m_particles.size() >= m_params.max_particles) {
        m_stats.skipped += stored;
        return;
    }
    m_particles.push_back({pos, vel, time, stored, bin});
}

void ChDustField::EmitFromWheels(const std::vector<WheelState>& wheels, double time, double dt) {
    const size_t nbins = m_params.bins.size();
    if (dt <= 0 || nbins == 0)
        return;
    if (m_carry.size() < wheels.size() * nbins)
        m_carry.resize(wheels.size() * nbins, 0.0);

    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::normal_distribution<double> normal(0.0, 1.0);
    const ChVector3d up(0, 0, 1);

    for (size_t w = 0; w < wheels.size(); ++w) {
        const WheelState& wheel = wheels[w];
        const ChVector3d a = wheel.axle.GetNormalized();
        // Down the wheel's plane: the ground normal with its component along the axle removed
        ChVector3d u = up - a * up.Dot(a);
        if (u.Length() < 1e-6)
            continue;
        u.Normalize();

        const ChVector3d& c = wheel.center;
        const double r = wheel.radius;
        const double ground = Ground(c.x(), c.y());
        const ChVector3d bottom = c - u * r;
        if (bottom.z() - ground > m_wheel.contact_tolerance)
            continue;

        const double rim_speed = std::abs(wheel.omega) * r;
        if (rim_speed < 1e-4)
            continue;

        // The rim at the bottom moves along -omega e relative to the center
        const ChVector3d e = a.Cross(u).GetNormalized();
        const ChVector3d d_bottom = e * (wheel.omega > 0 ? -1.0 : 1.0);
        const ChVector3d v_bottom = wheel.velocity - e * (wheel.omega * r);
        const double slip_speed = (v_bottom - u * v_bottom.Dot(u)).Length();
        const double mass_rate = m_wheel.bulk_density * wheel.width * m_wheel.loose_depth * (rim_speed + m_wheel.slip_gain * slip_speed);

        // Release angles, from the bottom of the wheel: from where the rim leaves the ground (the wheel sinks in by
        // r - (center height above ground)) up to the release limit
        const double phi_min = std::acos(std::clamp((c.z() - ground) / r, -1.0, 1.0));
        const double phi_max = std::max(phi_min, m_wheel.max_release_angle);

        for (size_t b = 0; b < nbins; ++b) {
            // Each bin gets super-particles in proportion to the light its grains block, so every super-particle
            // counts about as much in the rendered dust
            const double expected = m_wheel.particles_per_second * dt * m_bin_share[b];
            if (expected <= 0)
                continue;
            double& carry = m_carry[w * nbins + b];
            carry += expected;
            const int count = static_cast<int>(carry);
            carry -= count;
            const double mass = mass_rate * m_params.bins[b].mass_fraction * dt / expected;
            for (int n = 0; n < count; ++n) {
                const double phi = phi_min + (phi_max - phi_min) * unit(m_rng);
                const double lateral = (unit(m_rng) - 0.5) * wheel.width;
                const ChVector3d offset = (d_bottom * std::sin(phi) - u * std::cos(phi)) * r + a * lateral;
                const ChVector3d v_rim = wheel.velocity + a.Cross(offset) * wheel.omega;
                const double fraction = m_wheel.min_speed_fraction + (m_wheel.max_speed_fraction - m_wheel.min_speed_fraction) * unit(m_rng);
                const ChVector3d jitter(normal(m_rng), normal(m_rng), normal(m_rng));
                const ChVector3d vel = v_rim * fraction + jitter * (m_wheel.spread * rim_speed);
                // Released during the interval, from where the wheel has carried the release point by then
                const double t_release = dt * unit(m_rng);
                Emit(c + offset + wheel.velocity * t_release, vel, mass, static_cast<int>(b), time + t_release);
            }
        }
    }
}

void ChDustField::Update(double time) {
    m_stats.aloft = 0;
    size_t kept = 0;
    for (size_t n = 0; n < m_particles.size(); ++n) {
        const Particle& p = m_particles[n];
        const double t = time - p.time0;
        if (t > 0) {
            if (t > m_params.max_flight_time) {
                m_stats.expired += p.mass;
                continue;
            }
            // Landed: below the ground on the way down. A grain released below the surface of a sunken wheel is
            // still rising out of it.
            const ChVector3d pos = Position(p, time);
            const double vz = p.vel0.z() - m_params.gravity * t;
            if (vz < 0 && pos.z() < Ground(pos.x(), pos.y())) {
                m_stats.landed += p.mass;
                continue;
            }
        }
        m_stats.aloft += p.mass;
        m_particles[kept++] = p;
    }
    m_particles.resize(kept);
    m_stats.particles = kept;
}

void ChDustField::UpdateGrid(double time, const ChVector3d& center, const ChVector3d& sun_dir) {
    Grid& g = m_grid;
    const double h = m_spec.voxel;
    const int nx = std::max(2, m_spec.nx);
    const int ny = std::max(2, m_spec.ny);
    const int nz = std::max(2, m_spec.nz);
    const ChVector3d origin(std::floor((center.x() - 0.5 * nx * h) / h) * h, std::floor((center.y() - 0.5 * ny * h) / h) * h, std::floor((center.z() - m_spec.below) / h) * h);
    const bool resized = g.nx != nx || g.ny != ny || g.nz != nz || g.voxel != h;
    const size_t count = size_t(nx) * size_t(ny) * size_t(nz);
    g.origin = origin;
    g.voxel = h;
    g.nx = nx;
    g.ny = ny;
    g.nz = nz;
    g.extinction.assign(count, 0.f);
    g.sun_transmittance.resize(count);
    g.sun_visibility.resize(count);
    const ChVector3d L = sun_dir.GetNormalized();
    g.sun_dir = L;
    ++g.version;

    // Cloud-in-cell: each particle's extinction spread over the eight voxels around it
    const double inv_volume = 1.0 / (h * h * h);
    for (const Particle& p : m_particles) {
        if (p.time0 > time)
            continue;
        const ChVector3d q = (Position(p, time) - origin) / h - ChVector3d(0.5, 0.5, 0.5);
        const int i = static_cast<int>(std::floor(q.x()));
        const int j = static_cast<int>(std::floor(q.y()));
        const int k = static_cast<int>(std::floor(q.z()));
        if (i < -1 || j < -1 || k < -1 || i >= nx || j >= ny || k >= nz)
            continue;
        const double f[3] = {q.x() - i, q.y() - j, q.z() - k};
        const double sigma = p.mass * m_mass_extinction[p.bin] * inv_volume;
        for (int c = 0; c < 8; ++c) {
            const int ii = i + (c & 1), jj = j + ((c >> 1) & 1), kk = k + ((c >> 2) & 1);
            if (ii < 0 || jj < 0 || kk < 0 || ii >= nx || jj >= ny || kk >= nz)
                continue;
            const double w = ((c & 1) ? f[0] : 1 - f[0]) * (((c >> 1) & 1) ? f[1] : 1 - f[1]) * (((c >> 2) & 1) ? f[2] : 1 - f[2]);
            g.extinction[g.Index(ii, jj, kk)] += static_cast<float>(w * sigma);
        }
    }

    TraceTransmittance(g, L);
    if (resized || !m_vis_valid || (origin - m_vis_origin).Length() > 1e-9 || (L - m_vis_sun).Length() > 1e-6) {
        TraceVisibility(g, L, [this](double x, double y) { return Ground(x, y); }, m_spec.shadow_distance);
        m_vis_origin = origin;
        m_vis_sun = L;
        m_vis_valid = true;
    }
}

void ChDustField::TraceTransmittance(Grid& g, const ChVector3d& sun_dir) {
    const Sweep s = MakeSweep(sun_dir.GetNormalized(), g.voxel);
    const int n[3] = {g.nx, g.ny, g.nz};
    const int nb = n[s.b], nc = n[s.c];
    auto index = [&](int ia, int ib, int ic) {
        int ijk[3];
        ijk[s.a] = ia;
        ijk[s.b] = ib;
        ijk[s.c] = ic;
        return g.Index(ijk[0], ijk[1], ijk[2]);
    };
    // Extinction at a point of slice ia (voxel indices along b and c), bilinear, taking the outermost voxels' values
    // out to the faces of the grid and none past them
    auto extinction = [&](int ia, double qb, double qc) {
        if (qb < -0.5 || qc < -0.5 || qb > nb - 0.5 || qc > nc - 0.5)
            return 0.0;
        qb = std::clamp(qb, 0.0, nb - 1.0);
        qc = std::clamp(qc, 0.0, nc - 1.0);
        const int b0 = std::min(static_cast<int>(qb), nb - 2);
        const int c0 = std::min(static_cast<int>(qc), nc - 2);
        const double fb = qb - b0, fc = qc - c0;
        return (1 - fb) * (1 - fc) * g.extinction[index(ia, b0, c0)] + fb * (1 - fc) * g.extinction[index(ia, b0 + 1, c0)] + (1 - fb) * fc * g.extinction[index(ia, b0, c0 + 1)] +
               fb * fc * g.extinction[index(ia, b0 + 1, c0 + 1)];
    };

    // Rays from the Sun, one voxel apart, carried through the slices in the direction the light travels. Each
    // accumulates its own optical depth exactly, and each voxel interpolates once between the four rays around it,
    // so the sweep does not blur the shadows it traces. At slice t the ray stored at (m, l) passes through
    // (m + phase_b, l + phase_c), m from -1 to nb and l from -1 to nc; the phases advance by -off per slice, and
    // the storage shifts by whole voxels as they wrap.
    const int wb = nb + 2, wc = nc + 2;
    std::vector<double> tau(size_t(wb) * wc, 0.0), tau_next(tau.size());
    std::vector<double> last(tau.size(), 0.0), last_next(tau.size());  // extinction at each ray's previous sample
    const double shift_b = -s.off_b, shift_c = -s.off_c;               // ray displacement per slice (voxels)
    std::int64_t floor_b = 0, floor_c = 0;
    const bool par = util::parallelGrid(size_t(nb) * size_t(nc));

    for (int t = 0; t < n[s.a]; ++t) {
        const int ia = s.step > 0 ? n[s.a] - 1 - t : t;
        const std::int64_t fb_t = static_cast<std::int64_t>(std::floor(shift_b * t));
        const std::int64_t fc_t = static_cast<std::int64_t>(std::floor(shift_c * t));
        const double phase_b = shift_b * t - fb_t, phase_c = shift_c * t - fc_t;
        const int db = static_cast<int>(fb_t - floor_b), dc = static_cast<int>(fc_t - floor_c);
        floor_b = fb_t;
        floor_c = fc_t;

#pragma omp parallel for num_threads(util::parallelThreads()) schedule(static) if (par)
        for (int l = 0; l < wc; ++l) {
            for (int m = 0; m < wb; ++m) {
                const size_t r = size_t(m) + size_t(wb) * l;
                const double sigma = extinction(ia, m - 1 + phase_b, l - 1 + phase_c);
                // The ray now stored here was stored at (m - db, l - dc) at the previous slice. A ray entering the
                // grid through a side starts clear, as a ray entering at the first slice does.
                const int pm = m - db, pl = l - dc;
                if (t == 0 || pm < 0 || pl < 0 || pm >= wb || pl >= wc) {
                    tau_next[r] = 0.5 * sigma * s.length;
                } else {
                    const size_t pr = size_t(pm) + size_t(wb) * pl;
                    tau_next[r] = tau[pr] + 0.5 * (sigma + last[pr]) * s.length;
                }
                last_next[r] = sigma;
            }
        }
        tau.swap(tau_next);
        last.swap(last_next);

        // Each voxel of the slice lies between the rays stored at (ib - 1, ib) and (ic - 1, ic), offset by the
        // phases: indices ib and ib + 1 in the padded storage
#pragma omp parallel for num_threads(util::parallelThreads()) schedule(static) if (par)
        for (int ic = 0; ic < nc; ++ic) {
            for (int ib = 0; ib < nb; ++ib) {
                const size_t r00 = size_t(ib) + size_t(wb) * ic;  // ray (ib - 1, ic - 1)
                const double depth = phase_b * phase_c * tau[r00] + (1 - phase_b) * phase_c * tau[r00 + 1] + phase_b * (1 - phase_c) * tau[r00 + wb] +
                                     (1 - phase_b) * (1 - phase_c) * tau[r00 + wb + 1];
                g.sun_transmittance[index(ia, ib, ic)] = static_cast<float>(std::exp(-depth));
            }
        }
    }
}

void ChDustField::TraceVisibility(Grid& g, const ChVector3d& sun_dir, const HeightFunction& ground, double shadow_distance) {
    const ChVector3d L = sun_dir.GetNormalized();
    const int nx = g.nx, ny = g.ny, nz = g.nz;
    const double h = g.voxel;
    auto column_x = [&](int i) { return g.origin.x() + (i + 0.5) * h; };
    auto column_y = [&](int j) { return g.origin.y() + (j + 0.5) * h; };

    const double horizontal = std::sqrt(L.x() * L.x() + L.y() * L.y());
    if (L.z() <= 0) {
        // The Sun is below the horizon
        std::fill(g.sun_visibility.begin(), g.sun_visibility.end(), 0.f);
        return;
    }

    // Per column, the height below which the terrain hides the Sun
    std::vector<double> shade(size_t(nx) * size_t(ny));
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            shade[i + size_t(nx) * j] = ground(column_x(i), column_y(j));

    if (horizontal > 1e-9) {
        // Sweep the columns toward the Sun's azimuth: S = max(ground, S upstream - rise to it), with the upstream
        // column one column along the dominant horizontal axis away and interpolated along the other.
        const ChVector2d d(L.x() / horizontal, L.y() / horizontal);
        const double tan_elev = L.z() / horizontal;
        const bool along_x = std::abs(d.x()) >= std::abs(d.y());
        const int na = along_x ? nx : ny;
        const int nb = along_x ? ny : nx;
        const double da = along_x ? d.x() : d.y();
        const double db = along_x ? d.y() : d.x();
        const int step = da > 0 ? 1 : -1;
        const double run = h / std::abs(da);  // horizontal distance to the upstream column (m)
        const double off = db * run / h;      // its offset along the other axis (columns)
        const double march = std::max(h, 0.25);
        auto cell = [&](int ia, int ib) { return along_x ? size_t(ia) + size_t(nx) * ib : size_t(ib) + size_t(nx) * ia; };

        for (int k = 0; k < na; ++k) {
            const int ia = step > 0 ? na - 1 - k : k;
            const int qa = ia + step;
            for (int ib = 0; ib < nb; ++ib) {
                const size_t c = cell(ia, ib);
                const double qb = ib + off;
                double upstream;
                if (qa < 0 || qa >= na || qb < 0 || qb > nb - 1) {
                    // Past the grid: trace the ground toward the Sun
                    const double x = along_x ? column_x(ia) : column_x(ib);
                    const double y = along_x ? column_y(ib) : column_y(ia);
                    upstream = -1e30;
                    for (double t = run; t <= shadow_distance; t += march)
                        upstream = std::max(upstream, ground(x + d.x() * t, y + d.y() * t) - tan_elev * t);
                } else {
                    const int b0 = std::min(static_cast<int>(qb), nb - 2);
                    const double fb = qb - b0;
                    upstream = (1 - fb) * shade[cell(qa, b0)] + fb * shade[cell(qa, b0 + 1)] - tan_elev * run;
                }
                shade[c] = std::max(shade[c], upstream);
            }
        }
    }

    for (int k = 0; k < nz; ++k) {
        const double z = g.origin.z() + (k + 0.5) * h;
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const double v = (z - shade[i + size_t(nx) * j]) / h + 0.5;
                g.sun_visibility[g.Index(i, j, k)] = static_cast<float>(std::clamp(v, 0.0, 1.0));
            }
    }
}

bool ChDustField::WriteOpticalDepthImage(const std::string& filename) const {
    const Grid& g = m_grid;
    if (g.extinction.empty())
        return false;
    std::ofstream out(filename, std::ios::binary);
    if (!out)
        return false;
    out << "P5\n" << g.nx << " " << g.ny << "\n255\n";
    std::vector<unsigned char> row(g.nx);
    // North (+y) up
    for (int j = g.ny - 1; j >= 0; --j) {
        for (int i = 0; i < g.nx; ++i) {
            double tau = 0;
            for (int k = 0; k < g.nz; ++k)
                tau += g.extinction[g.Index(i, j, k)] * g.voxel;
            row[i] = static_cast<unsigned char>(std::lround(255.0 * (1.0 - std::exp(-tau))));
        }
        out.write(reinterpret_cast<const char*>(row.data()), row.size());
    }
    return static_cast<bool>(out);
}

}  // namespace planet
}  // namespace chrono
