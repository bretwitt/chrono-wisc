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

#include "chrono_planet/volume/ChSiteVolume.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <stdexcept>

#include "chrono/utils/ChConstants.h"

#include "chrono_planet/core/Parallel.h"
#include "chrono_planet/filters/ChDeformationFilter.h"

namespace chrono {
namespace planet {

namespace {

constexpr int B = ChSiteVolume::kBrick;

// Samples are kept to this many voxels from the ground
constexpr double kBandVoxels = 3.0;

int FloorDiv(int a, int b) {
    return a >= 0 ? a / b : -((-a + b - 1) / b);
}

size_t Local(int li, int lj, int lk) {
    return size_t(li) + size_t(B) * (size_t(lj) + size_t(B) * size_t(lk));
}

}  // namespace

ChSiteVolume::ChSiteVolume(const ChSiteRegion& region, const Params& params, HeightFunction ground)
    : m_params(params), m_h(params.voxel) {
    if (!(m_h > 0))
        throw std::invalid_argument("ChSiteVolume: voxel size must be positive");
    if (region.IsEmpty())
        throw std::invalid_argument("ChSiteVolume: empty region");
    if (!ground)
        throw std::invalid_argument("ChSiteVolume: no height function");
    m_band = static_cast<float>(kBandVoxels * m_h);

    // Lattice columns at whole voxels, covering the region with whole bricks
    const double x0 = std::floor(region.min_x / m_h) * m_h;
    const double y0 = std::floor(region.min_y / m_h) * m_h;
    for (int a = 0; a < 2; ++a) {
        const double span = a == 0 ? region.max_x - x0 : region.max_y - y0;
        m_nb[a] = std::max(1, static_cast<int>(std::ceil((span / m_h + 1) / B - 1e-9)));
        m_n[a] = m_nb[a] * B;
    }
    m_region = ChSiteRegion(x0, y0, x0 + (m_n[0] - 1) * m_h, y0 + (m_n[1] - 1) * m_h);
    m_edit_region = m_region.Inset(std::max(params.margin, 2 * m_h));

    // The initial ground, and the slope that turns its height offsets into distances
    m_ground.resize(size_t(m_n[0]) * m_n[1]);
    double gmin = 1e300, gmax = -1e300;
    for (int j = 0; j < m_n[1]; ++j)
        for (int i = 0; i < m_n[0]; ++i) {
            const double g = ground(x0 + i * m_h, y0 + j * m_h);
            m_ground[ColumnIndex(i, j)] = static_cast<float>(g);
            gmin = std::min(gmin, g);
            gmax = std::max(gmax, g);
        }
    m_slope.resize(m_ground.size());
    for (int j = 0; j < m_n[1]; ++j)
        for (int i = 0; i < m_n[0]; ++i) {
            const int il = std::max(0, i - 1), ir = std::min(m_n[0] - 1, i + 1);
            const int jl = std::max(0, j - 1), jr = std::min(m_n[1] - 1, j + 1);
            const double gx = (m_ground[ColumnIndex(ir, j)] - m_ground[ColumnIndex(il, j)]) / ((ir - il) * m_h);
            const double gy = (m_ground[ColumnIndex(i, jr)] - m_ground[ColumnIndex(i, jl)]) / ((jr - jl) * m_h);
            m_slope[ColumnIndex(i, j)] = static_cast<float>(1.0 / std::sqrt(1 + gx * gx + gy * gy));
        }

    const double z0 = std::floor((gmin - params.depth) / m_h) * m_h;
    const double z1 = std::ceil((gmax + params.height) / m_h) * m_h;
    m_nb[2] = std::max(1, static_cast<int>(std::ceil(((z1 - z0) / m_h + 1) / B - 1e-9)));
    m_n[2] = m_nb[2] * B;
    m_origin = ChVector3d(x0, y0, z0);

    // Bricks: those wholly above or below the band hold no samples
    const size_t num_bricks = size_t(m_nb[0]) * m_nb[1] * m_nb[2];
    m_table.assign(num_bricks, kAbove);
    m_brick_version.assign(num_bricks, 0);
    m_version = 1;
    const double h3 = m_h * m_h * m_h;
    double soil = 0;
    for (int bj = 0; bj < m_nb[1]; ++bj)
        for (int bi = 0; bi < m_nb[0]; ++bi) {
            double brick_gmin = 1e300, brick_gmax = -1e300, smin = 1;
            for (int lj = 0; lj < B; ++lj)
                for (int li = 0; li < B; ++li) {
                    const size_t c = ColumnIndex(bi * B + li, bj * B + lj);
                    brick_gmin = std::min(brick_gmin, double(m_ground[c]));
                    brick_gmax = std::max(brick_gmax, double(m_ground[c]));
                    smin = std::min(smin, double(m_slope[c]));
                }
            for (int bk = 0; bk < m_nb[2]; ++bk) {
                const size_t b = BrickIndex(bi, bj, bk);
                const double zlo = z0 + bk * B * m_h, zhi = zlo + (B - 1) * m_h;
                if ((zlo - brick_gmax) * smin >= m_band) {
                    m_table[b] = kAbove;
                    continue;
                }
                if ((zhi - brick_gmin) * smin <= -m_band) {
                    m_table[b] = kBelow;
                    soil += B * B * B;
                    continue;
                }
                Brick brick;
                bool above = true, below = true;
                for (int lk = 0; lk < B; ++lk)
                    for (int lj = 0; lj < B; ++lj)
                        for (int li = 0; li < B; ++li) {
                            const float d = InitialSample(bi * B + li, bj * B + lj, bk * B + lk);
                            brick[Local(li, lj, lk)] = d;
                            above = above && d >= m_band;
                            below = below && d <= -m_band;
                            soil += Fill(d);
                        }
                if (above || below) {
                    m_table[b] = above ? kAbove : kBelow;
                    continue;
                }
                m_table[b] = static_cast<std::int32_t>(m_pool.size());
                m_pool.push_back(brick);
                m_brick_version[b] = m_version;
            }
        }
    m_soil_volume = soil * h3;

    m_top.resize(m_ground.size());
    UpdateTops(0, 0, m_n[0] - 1, m_n[1] - 1);
    m_dirty_i0 = m_dirty_j0 = INT_MAX;
    m_dirty_i1 = m_dirty_j1 = INT_MIN;
}

// -----------------------------------------------------------------------------

float ChSiteVolume::InitialSample(int i, int j, int k) const {
    i = std::clamp(i, 0, m_n[0] - 1);
    j = std::clamp(j, 0, m_n[1] - 1);
    const size_t c = ColumnIndex(i, j);
    const double d = (m_origin.z() + k * m_h - m_ground[c]) * m_slope[c];
    return static_cast<float>(std::clamp(d, -double(m_band), double(m_band)));
}

float ChSiteVolume::Sample(int i, int j, int k) const {
    if (k < 0)
        return -m_band;
    if (k >= m_n[2])
        return m_band;
    i = std::clamp(i, 0, m_n[0] - 1);
    j = std::clamp(j, 0, m_n[1] - 1);
    const std::int32_t entry = m_table[BrickIndex(i / B, j / B, k / B)];
    if (entry == kAbove)
        return m_band;
    if (entry == kBelow)
        return -m_band;
    return m_pool[entry][Local(i % B, j % B, k % B)];
}

double ChSiteVolume::Fill(float d) const {
    return std::clamp(0.5 - d / m_h, 0.0, 1.0);
}

double ChSiteVolume::GetDistance(const ChVector3d& p) const {
    const ChVector3d u = (p - m_origin) / m_h;
    const int i = static_cast<int>(std::floor(u.x())), j = static_cast<int>(std::floor(u.y())), k = static_cast<int>(std::floor(u.z()));
    const double fx = u.x() - i, fy = u.y() - j, fz = u.z() - k;
    double d = 0;
    for (int c = 0; c < 8; ++c) {
        const int di = c & 1, dj = (c >> 1) & 1, dk = (c >> 2) & 1;
        const double w = (di ? fx : 1 - fx) * (dj ? fy : 1 - fy) * (dk ? fz : 1 - fz);
        if (w != 0)
            d += w * Sample(i + di, j + dj, k + dk);
    }
    return d;
}

ChVector3d ChSiteVolume::GetNormal(const ChVector3d& p) const {
    const double e = 0.5 * m_h;
    ChVector3d g(GetDistance(p + ChVector3d(e, 0, 0)) - GetDistance(p - ChVector3d(e, 0, 0)),
                 GetDistance(p + ChVector3d(0, e, 0)) - GetDistance(p - ChVector3d(0, e, 0)),
                 GetDistance(p + ChVector3d(0, 0, e)) - GetDistance(p - ChVector3d(0, 0, e)));
    const double len = g.Length();
    return len > 1e-12 ? g / len : ChVector3d(0, 0, 1);
}

double ChSiteVolume::ColumnTop(int i, int j) const {
    // Down the column from the top, skipping bricks of open space, to the first sample in the soil
    const int bi = i / B, bj = j / B, li = i % B, lj = j % B;
    for (int bk = m_nb[2] - 1; bk >= 0; --bk) {
        const std::int32_t entry = m_table[BrickIndex(bi, bj, bk)];
        if (entry == kAbove)
            continue;
        if (entry == kBelow)
            return m_origin.z() + (bk * B + B - 1) * m_h;
        const Brick& brick = m_pool[entry];
        for (int lk = B - 1; lk >= 0; --lk) {
            const float d = brick[Local(li, lj, lk)];
            if (d >= 0)
                continue;
            // The sample above is not in the soil, so the ground crosses between the two
            const int k = bk * B + lk;
            const float above = Sample(i, j, k + 1);
            return m_origin.z() + (k + d / (d - above)) * m_h;
        }
    }
    return m_origin.z();
}

void ChSiteVolume::UpdateTops(int i0, int j0, int i1, int j1) {
    for (int j = std::max(0, j0); j <= std::min(m_n[1] - 1, j1); ++j)
        for (int i = std::max(0, i0); i <= std::min(m_n[0] - 1, i1); ++i)
            m_top[ColumnIndex(i, j)] = static_cast<float>(ColumnTop(i, j));
}

namespace {
// Bilinear interpolation of a per-column field, clamped to the lattice
double Bilinear(const std::vector<float>& field, int nx, int ny, double u, double v) {
    u = std::clamp(u, 0.0, nx - 1.0);
    v = std::clamp(v, 0.0, ny - 1.0);
    const int i = std::min(static_cast<int>(u), nx - 2), j = std::min(static_cast<int>(v), ny - 2);
    const double fu = u - i, fv = v - j;
    auto at = [&](int a, int b) { return double(field[size_t(a) + size_t(nx) * size_t(b)]); };
    return (1 - fu) * (1 - fv) * at(i, j) + fu * (1 - fv) * at(i + 1, j) + (1 - fu) * fv * at(i, j + 1) + fu * fv * at(i + 1, j + 1);
}
}  // namespace

double ChSiteVolume::GetTopHeight(double x, double y) const {
    return Bilinear(m_top, m_n[0], m_n[1], (x - m_origin.x()) / m_h, (y - m_origin.y()) / m_h);
}

double ChSiteVolume::GetInitialHeight(double x, double y) const {
    return Bilinear(m_ground, m_n[0], m_n[1], (x - m_origin.x()) / m_h, (y - m_origin.y()) / m_h);
}

// -----------------------------------------------------------------------------

bool ChSiteVolume::EditRange(const ChAABB& box, int& i0, int& j0, int& k0, int& i1, int& j1, int& k1) const {
    if (box.IsInverted())
        return false;
    auto lo = [&](double v, double o) { return static_cast<int>(std::ceil((v - o) / m_h - 1e-9)); };
    auto hi = [&](double v, double o) { return static_cast<int>(std::floor((v - o) / m_h + 1e-9)); };
    i0 = std::max(lo(box.min.x(), m_origin.x()), lo(m_edit_region.min_x, m_origin.x()));
    i1 = std::min(hi(box.max.x(), m_origin.x()), hi(m_edit_region.max_x, m_origin.x()));
    j0 = std::max(lo(box.min.y(), m_origin.y()), lo(m_edit_region.min_y, m_origin.y()));
    j1 = std::min(hi(box.max.y(), m_origin.y()), hi(m_edit_region.max_y, m_origin.y()));
    k0 = std::max(lo(box.min.z(), m_origin.z()), 0);
    k1 = std::min(hi(box.max.z(), m_origin.z()), m_n[2] - 1);
    return i0 <= i1 && j0 <= j1 && k0 <= k1;
}

double ChSiteVolume::MeasureEdit(const ChAABB& box, const std::function<float(float, const ChVector3d&)>& op) const {
    int i0, j0, k0, i1, j1, k1;
    if (!EditRange(box, i0, j0, k0, i1, j1, k1))
        return 0;
    double fill_change = 0;
    for (int bk = k0 / B; bk <= k1 / B; ++bk)
        for (int bj = j0 / B; bj <= j1 / B; ++bj)
            for (int bi = i0 / B; bi <= i1 / B; ++bi) {
                const std::int32_t entry = m_table[BrickIndex(bi, bj, bk)];
                const Brick* brick = entry >= 0 ? &m_pool[entry] : nullptr;
                const float uniform = entry == kAbove ? m_band : -m_band;
                for (int k = std::max(k0, bk * B); k <= std::min(k1, bk * B + B - 1); ++k)
                    for (int j = std::max(j0, bj * B); j <= std::min(j1, bj * B + B - 1); ++j)
                        for (int i = std::max(i0, bi * B); i <= std::min(i1, bi * B + B - 1); ++i) {
                            const float old = brick ? (*brick)[Local(i - bi * B, j - bj * B, k - bk * B)] : uniform;
                            const float now = std::clamp(op(old, PointPosition(i, j, k)), -m_band, m_band);
                            if (now != old)
                                fill_change += Fill(now) - Fill(old);
                        }
            }
    return fill_change * m_h * m_h * m_h;
}

double ChSiteVolume::Edit(const ChAABB& box, const std::function<float(float, const ChVector3d&)>& op) {
    int i0, j0, k0, i1, j1, k1;
    if (!EditRange(box, i0, j0, k0, i1, j1, k1))
        return 0;

    const std::uint64_t version = m_version + 1;
    double fill_change = 0;
    int ci0 = INT_MAX, cj0 = INT_MAX, ck0 = INT_MAX, ci1 = INT_MIN, cj1 = INT_MIN, ck1 = INT_MIN;
    Brick scratch;
    for (int bk = k0 / B; bk <= k1 / B; ++bk)
        for (int bj = j0 / B; bj <= j1 / B; ++bj)
            for (int bi = i0 / B; bi <= i1 / B; ++bi) {
                const size_t b = BrickIndex(bi, bj, bk);
                const std::int32_t entry = m_table[b];
                if (entry >= 0)
                    scratch = m_pool[entry];
                else
                    scratch.fill(entry == kAbove ? m_band : -m_band);

                bool changed = false;
                for (int k = std::max(k0, bk * B); k <= std::min(k1, bk * B + B - 1); ++k)
                    for (int j = std::max(j0, bj * B); j <= std::min(j1, bj * B + B - 1); ++j)
                        for (int i = std::max(i0, bi * B); i <= std::min(i1, bi * B + B - 1); ++i) {
                            float& d = scratch[Local(i - bi * B, j - bj * B, k - bk * B)];
                            const float old = d;
                            const float now = std::clamp(op(old, PointPosition(i, j, k)), -m_band, m_band);
                            if (now == old)
                                continue;
                            d = now;
                            fill_change += Fill(now) - Fill(old);
                            changed = true;
                            ci0 = std::min(ci0, i), ci1 = std::max(ci1, i);
                            cj0 = std::min(cj0, j), cj1 = std::max(cj1, j);
                            ck0 = std::min(ck0, k), ck1 = std::max(ck1, k);
                        }
                if (!changed)
                    continue;

                // Store the brick, or drop its samples if they are all open space or all soil
                bool above = true, below = true;
                for (float d : scratch) {
                    above = above && d >= m_band;
                    below = below && d <= -m_band;
                }
                if (above || below) {
                    if (entry >= 0)
                        m_free.push_back(entry);
                    m_table[b] = above ? kAbove : kBelow;
                } else if (entry >= 0) {
                    m_pool[entry] = scratch;
                } else if (!m_free.empty()) {
                    m_table[b] = m_free.back();
                    m_free.pop_back();
                    m_pool[m_table[b]] = scratch;
                } else {
                    m_table[b] = static_cast<std::int32_t>(m_pool.size());
                    m_pool.push_back(scratch);
                }
            }
    if (ci0 > ci1)
        return 0;

    // A changed point moves the ground in the cells around it, which the bricks meshing those cells hold
    m_version = version;
    for (int bk = std::max(0, FloorDiv(ck0 - 1, B)); bk <= std::min(m_nb[2] - 1, FloorDiv(ck1 + 1, B)); ++bk)
        for (int bj = std::max(0, FloorDiv(cj0 - 1, B)); bj <= std::min(m_nb[1] - 1, FloorDiv(cj1 + 1, B)); ++bj)
            for (int bi = std::max(0, FloorDiv(ci0 - 1, B)); bi <= std::min(m_nb[0] - 1, FloorDiv(ci1 + 1, B)); ++bi)
                m_brick_version[BrickIndex(bi, bj, bk)] = version;

    UpdateTops(ci0, cj0, ci1, cj1);
    m_dirty_i0 = std::min(m_dirty_i0, ci0), m_dirty_i1 = std::max(m_dirty_i1, ci1);
    m_dirty_j0 = std::min(m_dirty_j0, cj0), m_dirty_j1 = std::max(m_dirty_j1, cj1);

    const double change = fill_change * m_h * m_h * m_h;
    m_soil_volume += change;
    return change;
}

double ChSiteVolume::Subtract(const ChSdfShape& shape, const ChFrame<>& frame) {
    const ChAABB box = shape.GetBoundingBox().Transform(frame);
    return -Edit(ChAABB(box.min - ChVector3d(m_band), box.max + ChVector3d(m_band)), [&](float d, const ChVector3d& p) {
        return std::max(d, static_cast<float>(-shape.Distance(frame.TransformPointParentToLocal(p))));
    });
}

double ChSiteVolume::Add(const ChSdfShape& shape, const ChFrame<>& frame) {
    const ChAABB box = shape.GetBoundingBox().Transform(frame);
    return Edit(ChAABB(box.min - ChVector3d(m_band), box.max + ChVector3d(m_band)), [&](float d, const ChVector3d& p) {
        return std::min(d, static_cast<float>(shape.Distance(frame.TransformPointParentToLocal(p))));
    });
}

double ChSiteVolume::SubtractSwept(const ChSdfShape& shape, const ChFrame<>& from, const ChFrame<>& to) {
    // Placements between the two, no point of the shape's box moving more than half a voxel from one to the next
    const ChAABB local = shape.GetBoundingBox();
    double travel = 0;
    for (int c = 0; c < 8; ++c) {
        const ChVector3d corner((c & 1) ? local.max.x() : local.min.x(), (c & 2) ? local.max.y() : local.min.y(), (c & 4) ? local.max.z() : local.min.z());
        travel = std::max(travel, (to.TransformPointLocalToParent(corner) - from.TransformPointLocalToParent(corner)).Length());
    }
    const int steps = std::max(1, static_cast<int>(std::ceil(travel / (0.5 * m_h))));
    ChQuaternion<> q0 = from.GetRot(), q1 = to.GetRot();
    if (q0.Dot(q1) < 0)
        q1 = -q1;
    std::vector<ChFrame<>> frames;
    frames.reserve(steps + 1);
    ChAABB box;
    for (int s = 0; s <= steps; ++s) {
        const double t = double(s) / steps;
        ChQuaternion<> q = q0 * (1 - t) + q1 * t;
        q.Normalize();
        frames.emplace_back(from.GetPos() * (1 - t) + to.GetPos() * t, q);
        box += local.Transform(frames.back());
    }
    return -Edit(ChAABB(box.min - ChVector3d(m_band), box.max + ChVector3d(m_band)), [&](float d, const ChVector3d& p) {
        double inside = -1e30;
        for (const auto& f : frames)
            inside = std::max(inside, -shape.Distance(f.TransformPointParentToLocal(p)));
        return std::max(d, static_cast<float>(inside));
    });
}

double ChSiteVolume::Deposit(double x, double y, double volume, double repose_angle) {
    if (!(volume > 0))
        return 0;
    const double angle = std::clamp(repose_angle, 0.05, 1.5);
    const double slope = std::tan(angle), cos_angle = std::cos(angle);

    // The lowest top the pile could reach down to, within twice the reach of the same pile on flat ground, which
    // bounds how wide it spreads and how deep it fills
    auto lo = [&](double v, double o) { return static_cast<int>(std::ceil((v - o) / m_h - 1e-9)); };
    auto hi = [&](double v, double o) { return static_cast<int>(std::floor((v - o) / m_h + 1e-9)); };
    const double base = GetTopHeight(x, y);
    const double flat_reach = std::cbrt(3 * volume / (CH_PI * slope));
    const double far = 2 * flat_reach + 4 * m_h;
    const int ci0 = std::max(0, lo(x - far, m_origin.x())), ci1 = std::min(m_n[0] - 1, hi(x + far, m_origin.x()));
    const int cj0 = std::max(0, lo(y - far, m_origin.y())), cj1 = std::min(m_n[1] - 1, hi(y + far, m_origin.y()));
    double lowest = base;
    for (int j = cj0; j <= cj1; ++j)
        for (int i = ci0; i <= ci1; ++i)
            lowest = std::min(lowest, double(m_top[ColumnIndex(i, j)]));
    auto reach = [&](double a) { return std::min(far, std::max(0.0, (a - lowest) / slope)); };

    // A first apex from the pile's volume over the column tops, cheap to evaluate
    auto over_columns = [&](double a) {
        const double r = reach(a);
        double v = 0;
        for (int j = std::max(cj0, lo(y - r, m_origin.y())); j <= std::min(cj1, hi(y + r, m_origin.y())); ++j)
            for (int i = std::max(ci0, lo(x - r, m_origin.x())); i <= std::min(ci1, hi(x + r, m_origin.x())); ++i) {
                const ChVector3d p = PointPosition(i, j, 0);
                const double dx = p.x() - x, dy = p.y() - y;
                v += std::max(0.0, a - std::sqrt(dx * dx + dy * dy) * slope - m_top[ColumnIndex(i, j)]);
            }
        return v * m_h * m_h;
    };
    double a_lo = lowest, a_hi = base + slope * flat_reach + m_h;
    for (int n = 0; n < 20 && over_columns(a_hi) < volume && a_hi < GetMaxZ(); ++n)
        a_hi = base + 2 * (a_hi - base);
    for (int n = 0; n < 40; ++n) {
        const double a = 0.5 * (a_lo + a_hi);
        (over_columns(a) < volume ? a_lo : a_hi) = a;
    }
    const double estimate = a_hi;

    // The union with a cone of apex a, over the box the cone can change: its reach, and from the band below the
    // lowest top up to the apex
    auto box_of = [&](double a) {
        const double r = reach(a) + m_h;
        return ChAABB(ChVector3d(x - r, y - r, lowest - m_band - m_h), ChVector3d(x + r, y + r, a + m_band));
    };
    auto cone = [&](double a) {
        return [=](float d, const ChVector3d& p) {
            const double dx = p.x() - x, dy = p.y() - y;
            const double radius = std::sqrt(dx * dx + dy * dy);
            return std::min(d, static_cast<float>((p.z() - (a - radius * slope)) * cos_angle));
        };
    };
    // The lattice counts the soil a little differently from the columns, mostly the half voxel under the old
    // ground that the pile covers: refine the apex by the volume the edit will measure
    auto added = [&](double a) { return MeasureEdit(box_of(a), cone(a)); };
    a_lo = estimate - 2 * m_h;
    for (int n = 0; n < 20 && added(a_lo) > volume && a_lo > lowest - m_band; ++n)
        a_lo -= 2 * m_h;
    a_hi = estimate + m_h;
    for (int n = 0; n < 20 && added(a_hi) < volume && a_hi < GetMaxZ(); ++n)
        a_hi += 2 * (a_hi - estimate);
    for (int n = 0; n < 12; ++n) {
        const double a = 0.5 * (a_lo + a_hi);
        (added(a) < volume ? a_lo : a_hi) = a;
    }
    const double apex = std::abs(added(a_lo) - volume) < std::abs(added(a_hi) - volume) ? a_lo : a_hi;
    return Edit(box_of(apex), cone(apex));
}

double ChSiteVolume::Assign(const ChAABB& box, const std::function<float(const ChVector3d& p)>& distance) {
    return Edit(box, [&](float, const ChVector3d& p) { return distance(p); });
}

// -----------------------------------------------------------------------------

void ChSiteVolume::GetChangedBricks(std::uint64_t since, std::vector<ChVector3i>& bricks) const {
    bricks.clear();
    for (int bk = 0; bk < m_nb[2]; ++bk)
        for (int bj = 0; bj < m_nb[1]; ++bj)
            for (int bi = 0; bi < m_nb[0]; ++bi)
                if (m_brick_version[BrickIndex(bi, bj, bk)] > since)
                    bricks.push_back(ChVector3i(bi, bj, bk));
}

bool ChSiteVolume::MeshBrick(const ChVector3i& brick, ChTriangleMeshConnected& mesh, double disturbed) const {
    mesh.Clear();
    const int bi = brick.x(), bj = brick.y(), bk = brick.z();
    if (bi < 0 || bj < 0 || bk < 0 || bi >= m_nb[0] || bj >= m_nb[1] || bk >= m_nb[2])
        return false;
    if (m_table[BrickIndex(bi, bj, bk)] < 0)
        return false;

    // Samples at points [8b - 1, 8b + 8] along each axis, and whether they differ from the initial ground
    constexpr int S = B + 2;
    const int pi0 = bi * B - 1, pj0 = bj * B - 1, pk0 = bk * B - 1;
    std::array<float, S * S * S> d;
    std::array<bool, S * S * S> moved;
    auto sidx = [](int a, int b, int c) { return size_t(a) + size_t(S) * (size_t(b) + size_t(S) * size_t(c)); };
    for (int c = 0; c < S; ++c)
        for (int b = 0; b < S; ++b)
            for (int a = 0; a < S; ++a) {
                const float v = Sample(pi0 + a, pj0 + b, pk0 + c);
                d[sidx(a, b, c)] = v;
                moved[sidx(a, b, c)] = std::abs(v - InitialSample(pi0 + a, pj0 + b, pk0 + c)) > disturbed;
            }

    // Cells [8b - 1, 8b + 7], local index 0..8, each spanning samples c..c+1 of the cache
    constexpr int C = B + 1;
    auto cidx = [](int a, int b, int c) { return size_t(a) + size_t(C) * (size_t(b) + size_t(C) * size_t(c)); };
    auto cell_valid = [&](int a, int b, int c) {
        return pi0 + a >= 0 && pj0 + b >= 0 && pk0 + c >= 0 && pi0 + a + 1 < m_n[0] && pj0 + b + 1 < m_n[1] && pk0 + c + 1 < m_n[2];
    };
    std::array<int, C * C * C> vertex;
    vertex.fill(-2);  // -2: not computed, -1: no vertex

    auto& vertices = mesh.GetCoordsVertices();
    auto& normals = mesh.GetCoordsNormals();
    std::vector<bool> vertex_moved;

    // Trilinear distance in the cache at fractional sample coordinates, for normals
    auto cached = [&](double a, double b, double c) {
        a = std::clamp(a, 0.0, S - 1.0001), b = std::clamp(b, 0.0, S - 1.0001), c = std::clamp(c, 0.0, S - 1.0001);
        const int ia = static_cast<int>(a), ib = static_cast<int>(b), ic = static_cast<int>(c);
        const double fa = a - ia, fb = b - ib, fc = c - ic;
        double v = 0;
        for (int n = 0; n < 8; ++n) {
            const int da = n & 1, db = (n >> 1) & 1, dc = (n >> 2) & 1;
            v += (da ? fa : 1 - fa) * (db ? fb : 1 - fb) * (dc ? fc : 1 - fc) * d[sidx(ia + da, ib + db, ic + dc)];
        }
        return v;
    };

    // Vertex of a cell: the mean of the ground's crossings of its edges
    auto cell_vertex = [&](int a, int b, int c) -> int {
        int& slot = vertex[cidx(a, b, c)];
        if (slot != -2)
            return slot;
        slot = -1;
        if (!cell_valid(a, b, c))
            return slot;
        float corner[8];
        bool any_in = false, any_out = false, any_moved = false;
        for (int n = 0; n < 8; ++n) {
            const size_t s = sidx(a + (n & 1), b + ((n >> 1) & 1), c + ((n >> 2) & 1));
            corner[n] = d[s];
            (corner[n] < 0 ? any_in : any_out) = true;
            any_moved = any_moved || moved[s];
        }
        if (!any_in || !any_out)
            return slot;
        static const int edges[12][2] = {{0, 1}, {2, 3}, {4, 5}, {6, 7}, {0, 2}, {1, 3}, {4, 6}, {5, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
        ChVector3d sum(0, 0, 0);
        int count = 0;
        for (const auto& e : edges) {
            const float d0 = corner[e[0]], d1 = corner[e[1]];
            if ((d0 < 0) == (d1 < 0))
                continue;
            const double t = d0 / (d0 - d1);
            const ChVector3d p0((e[0] & 1), (e[0] >> 1) & 1, (e[0] >> 2) & 1);
            const ChVector3d p1((e[1] & 1), (e[1] >> 1) & 1, (e[1] >> 2) & 1);
            sum += p0 + (p1 - p0) * t;
            ++count;
        }
        ChVector3d f = sum / count;  // within the cell, in voxels
        // The mean of the crossings lies off the surface where it bends, which shows as terraces on slopes: step it
        // onto the zero of the trilinear field along its gradient, staying in the cell so faces cannot fold. The field
        // is read from the whole volume at the cell's lattice position, so bricks sharing the cell place it alike.
        const ChVector3d cell_origin = PointPosition(pi0 + a, pj0 + b, pk0 + c);
        const double e = 0.5 * m_h;
        for (int iter = 0; iter < 2; ++iter) {
            const ChVector3d p = cell_origin + f * m_h;
            const ChVector3d grad(GetDistance(p + ChVector3d(e, 0, 0)) - GetDistance(p - ChVector3d(e, 0, 0)),
                                  GetDistance(p + ChVector3d(0, e, 0)) - GetDistance(p - ChVector3d(0, e, 0)),
                                  GetDistance(p + ChVector3d(0, 0, e)) - GetDistance(p - ChVector3d(0, 0, e)));
            const double g2 = grad.Length2();
            if (g2 < 1e-20)
                break;
            f -= grad * (GetDistance(p) * m_h / g2);  // grad is per voxel: a step of phi / |grad| voxels
            for (int k = 0; k < 3; ++k)
                f[k] = std::clamp(f[k], 0.0, 1.0);
        }
        const double sa = a + f.x(), sb = b + f.y(), sc = c + f.z();
        ChVector3d g(cached(sa + 0.5, sb, sc) - cached(sa - 0.5, sb, sc), cached(sa, sb + 0.5, sc) - cached(sa, sb - 0.5, sc),
                     cached(sa, sb, sc + 0.5) - cached(sa, sb, sc - 0.5));
        const double len = g.Length();
        slot = static_cast<int>(vertices.size());
        vertices.push_back(PointPosition(pi0 + a, pj0 + b, pk0 + c) + f * m_h);
        normals.push_back(len > 1e-12 ? g / len : ChVector3d(0, 0, 1));
        vertex_moved.push_back(any_moved);
        return slot;
    };

    auto& faces = mesh.GetIndicesVertices();
    auto& face_normals = mesh.GetIndicesNormals();
    auto& face_materials = mesh.GetIndicesMaterials();
    auto add_triangle = [&](int v0, int v1, int v2) {
        faces.push_back(ChVector3i(v0, v1, v2));
        face_normals.push_back(ChVector3i(v0, v1, v2));
        face_materials.push_back(vertex_moved[v0] || vertex_moved[v1] || vertex_moved[v2] ? 1 : 0);
    };

    // A quad around each edge the ground crosses, for the edges starting at this brick's points
    for (int lk = 0; lk < B; ++lk)
        for (int lj = 0; lj < B; ++lj)
            for (int li = 0; li < B; ++li) {
                const int p[3] = {li + 1, lj + 1, lk + 1};  // in the cache
                const bool inside = d[sidx(p[0], p[1], p[2])] < 0;
                for (int axis = 0; axis < 3; ++axis) {
                    int q[3] = {p[0], p[1], p[2]};
                    ++q[axis];
                    if ((axis == 0 ? pi0 : axis == 1 ? pj0 : pk0) + q[axis] >= m_n[axis])
                        continue;
                    if ((d[sidx(q[0], q[1], q[2])] < 0) == inside)
                        continue;
                    // The four cells around the edge, cell index = sample index of its low corner
                    const int u = (axis + 1) % 3, v = (axis + 2) % 3;
                    int cells[4][3];
                    for (int n = 0; n < 4; ++n) {
                        for (int a = 0; a < 3; ++a)
                            cells[n][a] = p[a];
                        if (n == 1 || n == 2)
                            --cells[n][u];
                        if (n == 2 || n == 3)
                            --cells[n][v];
                    }
                    int ids[4];
                    bool ok = true;
                    for (int n = 0; n < 4 && ok; ++n) {
                        ids[n] = cell_vertex(cells[n][0], cells[n][1], cells[n][2]);
                        ok = ids[n] >= 0;
                    }
                    if (!ok)
                        continue;
                    // Counter-clockwise seen from outside the soil: the soil is on the p side of the edge if inside
                    if (!inside)
                        std::swap(ids[1], ids[3]);
                    if ((vertices[ids[0]] - vertices[ids[2]]).Length2() <= (vertices[ids[1]] - vertices[ids[3]]).Length2()) {
                        add_triangle(ids[0], ids[1], ids[2]);
                        add_triangle(ids[0], ids[2], ids[3]);
                    } else {
                        add_triangle(ids[0], ids[1], ids[3]);
                        add_triangle(ids[1], ids[2], ids[3]);
                    }
                }
            }
    return !faces.empty();
}

// -----------------------------------------------------------------------------

size_t ChSiteVolume::PublishTopSurface(ChDeformationFilter& filter, const ChSiteFrame& site) {
    if (std::abs(filter.GetSpacing() - m_h) > 1e-9 * m_h)
        throw std::invalid_argument("ChSiteVolume::PublishTopSurface: the filter's node spacing must equal the voxel size");
    if (m_dirty_i0 > m_dirty_i1)
        return 0;
    const int oi = static_cast<int>(std::lround(m_origin.x() / m_h));
    const int oj = static_cast<int>(std::lround(m_origin.y() / m_h));
    std::vector<ChDeformationFilter::Node> nodes;
    nodes.reserve(size_t(m_dirty_i1 - m_dirty_i0 + 1) * (m_dirty_j1 - m_dirty_j0 + 1));
    for (int j = m_dirty_j0; j <= m_dirty_j1; ++j)
        for (int i = m_dirty_i0; i <= m_dirty_i1; ++i) {
            const size_t c = ColumnIndex(i, j);
            nodes.push_back({ChVector2i(oi + i, oj + j), double(m_top[c]) - m_ground[c], double(m_top[c]) + site.GetOriginElevation()});
        }
    m_dirty_i0 = m_dirty_j0 = INT_MAX;
    m_dirty_i1 = m_dirty_j1 = INT_MIN;
    return filter.SetNodes(nodes, 1e-4);
}

}  // namespace planet
}  // namespace chrono
