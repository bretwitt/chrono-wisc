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

#include "chrono_planet/volume/ChSoilParticles.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

#include "chrono_planet/volume/ChSdfMesher.h"

namespace chrono {
namespace planet {

namespace {
// Radius of a clod, in spacings: from this to this plus kRadiusSpread, so clods of a pour differ in size
constexpr double kRadius = 0.45;
constexpr double kRadiusSpread = 0.2;

// Value noise in [-1, 1] over the plane, of unit wavelength, for the lumps of a load
double Hash(int i, int j) {
    std::uint32_t h = static_cast<std::uint32_t>(i) * 374761393u + static_cast<std::uint32_t>(j) * 668265263u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return ((h ^ (h >> 16)) & 0xffffff) / double(0xffffff) * 2 - 1;
}
double Noise(double x, double y) {
    const int i = static_cast<int>(std::floor(x)), j = static_cast<int>(std::floor(y));
    const double u = x - i, v = y - j;
    const double su = u * u * (3 - 2 * u), sv = v * v * (3 - 2 * v);
    return (1 - su) * (1 - sv) * Hash(i, j) + su * (1 - sv) * Hash(i + 1, j) + (1 - su) * sv * Hash(i, j + 1) + su * sv * Hash(i + 1, j + 1);
}
// Mesh resolution, in spacings
constexpr double kMeshStep = 0.5;
}  // namespace

ChSoilParticles::ChSoilParticles(ChSystem* sys, std::shared_ptr<ChSiteVolume> volume, const Params& params)
    : m_system(sys),
      m_volume(std::move(volume)),
      m_params(params),
      m_grid(kMeshStep * params.spacing, static_cast<float>((kRadius + params.blend + 2 * kMeshStep) * params.spacing)) {
    if (!m_system || !m_volume)
        throw std::invalid_argument("ChSoilParticles: null system or volume");
    if (!(m_params.spacing > 0))
        throw std::invalid_argument("ChSoilParticles: spacing must be positive");
    m_body = chrono_types::make_shared<ChBody>();
    m_body->SetFixed(true);
    m_body->EnableCollision(false);
    m_system->AddBody(m_body);
    m_material = chrono_types::make_shared<ChVisualMaterial>();
    m_material->SetDiffuseColor(ChColor(0.33f, 0.32f, 0.30f));
    m_material->SetRoughness(0.95f);
    m_material->SetMetallic(0.0f);
    m_mesh = chrono_types::make_shared<ChTriangleMeshConnected>();
}

ChSoilParticles::ChSoilParticles(ChSystem* sys, std::shared_ptr<ChSiteVolume> volume) : ChSoilParticles(sys, std::move(volume), Params()) {}

ChSoilParticles::~ChSoilParticles() {}

// -----------------------------------------------------------------------------

int ChSoilParticles::AddCarrier(std::shared_ptr<ChBody> body, const ChFrame<>& frame, const ChAABB& cavity) {
    if (!body || cavity.IsInverted())
        throw std::invalid_argument("ChSoilParticles::AddCarrier: null body or empty cavity");
    Carrier c;
    c.body = std::move(body);
    c.frame = frame;
    c.cavity = cavity;
    const ChVector3d size = cavity.Size();
    c.nx = std::max(1, static_cast<int>(size.x() / m_params.spacing));
    c.ny = std::max(1, static_cast<int>(size.y() / m_params.spacing));
    m_carriers.push_back(std::move(c));
    return static_cast<int>(m_carriers.size()) - 1;
}

ChVector3d ChSoilParticles::SlotPosition(const Carrier& c, int slot) const {
    const int a = slot % c.nx, b = (slot / c.nx) % c.ny, layer = slot / (c.nx * c.ny);
    // Center the lattice across the cavity; layers stack up from its floor, above its rim once it is full
    const ChVector3d size = c.cavity.Size();
    const double s = m_params.spacing;
    const double x0 = c.cavity.min.x() + 0.5 * (size.x() - c.nx * s), y0 = c.cavity.min.y() + 0.5 * (size.y() - c.ny * s);
    return ChVector3d(x0 + (a + 0.5) * s, y0 + (b + 0.5) * s, c.cavity.min.z() + (layer + 0.5) * s);
}

void ChSoilParticles::SetCarried(int carrier, double volume) {
    Carrier& c = m_carriers.at(carrier);
    const size_t want = static_cast<size_t>(std::max(0.0, std::floor(volume / GetParticleVolume() + 1e-9)));
    while (c.slots.size() > want)
        c.slots.pop_back();
    if (c.slots.size() == want)
        return;
    // New particles take the lowest free slots
    std::vector<bool> used;
    for (int s : c.slots) {
        if (s >= static_cast<int>(used.size()))
            used.resize(s + 1, false);
        used[s] = true;
    }
    int next = 0;
    while (c.slots.size() < want) {
        while (next < static_cast<int>(used.size()) && used[next])
            ++next;
        c.slots.push_back(next++);
    }
}

double ChSoilParticles::Release(int carrier, double volume) {
    Carrier& c = m_carriers.at(carrier);
    const size_t count = std::min(c.slots.size(), static_cast<size_t>(std::floor(volume / GetParticleVolume() + 1e-9)));
    if (count == 0)
        return 0;
    // The lowest particles go first: over the lip of a tipped bucket, and out of the bottom of an opened one
    const ChFrame<> frame = CarrierFrame(c);
    std::vector<std::pair<double, size_t>> order;
    for (size_t n = 0; n < c.slots.size(); ++n)
        order.push_back({frame.TransformPointLocalToParent(SlotPosition(c, c.slots[n])).z(), n});
    std::partial_sort(order.begin(), order.begin() + count, order.end());
    std::vector<bool> gone(c.slots.size(), false);
    for (size_t n = 0; n < count; ++n) {
        const size_t idx = order[n].second;
        std::uniform_real_distribution<double> unit(-1, 1);
        std::normal_distribution<double> normal(0, m_params.spread);
        const ChVector3d jitter(unit(m_rng), unit(m_rng), unit(m_rng));
        const ChVector3d pos = frame.TransformPointLocalToParent(SlotPosition(c, c.slots[idx]) + jitter * (0.4 * m_params.spacing));
        const ChVector3d vel = c.body->PointSpeedLocalToParent(c.body->TransformPointParentToLocal(pos)) + ChVector3d(normal(m_rng), normal(m_rng), normal(m_rng));
        m_free.push_back({pos, vel, false, static_cast<float>(kRadius + kRadiusSpread * std::uniform_real_distribution<double>(0, 1)(m_rng))});
        gone[idx] = true;
    }
    std::vector<int> kept;
    for (size_t n = 0; n < c.slots.size(); ++n)
        if (!gone[n])
            kept.push_back(c.slots[n]);
    c.slots.swap(kept);
    return count * GetParticleVolume();
}

double ChSoilParticles::GetCarriedVolume(int carrier) const {
    return m_carriers.at(carrier).slots.size() * GetParticleVolume();
}

size_t ChSoilParticles::GetNumCarried() const {
    size_t n = 0;
    for (const auto& c : m_carriers)
        n += c.slots.size();
    return n;
}

size_t ChSoilParticles::GetNumFalling() const {
    return std::count_if(m_free.begin(), m_free.end(), [](const Particle& p) { return !p.resting; });
}

size_t ChSoilParticles::GetNumResting() const {
    return std::count_if(m_free.begin(), m_free.end(), [](const Particle& p) { return p.resting; });
}

// -----------------------------------------------------------------------------

void ChSoilParticles::Advance(double dt) {
    if (!(dt > 0))
        return;
    const double r = kRadius * m_params.spacing;
    for (auto& p : m_free) {
        if (p.resting)
            continue;
        p.vel.z() -= m_params.gravity * dt;
        p.pos += p.vel * dt;
        // Resting once the particle's sphere touches the ground, set down on it
        const double d = m_volume->GetDistance(p.pos);
        if (d < r) {
            p.pos += m_volume->GetNormal(p.pos) * (r - d);
            p.vel = VNULL;
            p.resting = true;
        }
    }
    m_since_merge += dt;
    if (m_since_merge >= m_params.merge_period)
        Merge();
}

void ChSoilParticles::Merge() {
    m_since_merge = 0;
    // Resting particles by cell, each cell's poured as one pile of their volume at their mean position
    struct Cell {
        double x = 0, y = 0;
        int count = 0;
    };
    std::map<std::pair<long, long>, Cell> cells;
    const double w = m_params.merge_cell;
    std::vector<Particle> falling;
    for (const auto& p : m_free) {
        if (!p.resting) {
            falling.push_back(p);
            continue;
        }
        Cell& cell = cells[{static_cast<long>(std::floor(p.pos.x() / w)), static_cast<long>(std::floor(p.pos.y() / w))}];
        cell.x += p.pos.x();
        cell.y += p.pos.y();
        ++cell.count;
    }
    if (cells.empty())
        return;
    for (const auto& [key, cell] : cells) {
        const double volume = cell.count * GetParticleVolume();
        const double added = m_volume->Deposit(cell.x / cell.count, cell.y / cell.count, volume, m_params.repose_angle);
        m_merged += added;
        m_lost += volume - added;
    }
    m_free.swap(falling);
}

// -----------------------------------------------------------------------------

void ChSoilParticles::UpdateMesh() {
    // Everything out of the ground, splatted where it is now; bricks whose samples did not change keep their mesh
    const double s = m_params.spacing;
    const double blend = std::max(0.0, m_params.blend * s);
    m_grid.Begin();

    // Each carrier's load, as the soil filling its cavity: level where it holds the load, a cone at the angle of
    // repose above the rim for what does not fit, and lumps on top, fixed in the carrier's frame
    const double slope = std::tan(m_params.repose_angle);
    for (const auto& c : m_carriers) {
        if (c.slots.empty())
            continue;
        const ChVector3d size = c.cavity.Size();
        const ChVector3d center = c.cavity.Center();
        const double area = size.x() * size.y();
        const double load = c.slots.size() * GetParticleVolume();
        const double level = std::min(size.z(), load / area);
        const double over = std::max(0.0, load - area * size.z());
        const double cone = over > 0 ? std::cbrt(3 * over * slope * slope / CH_PI) : 0;
        const double lumps = m_params.lumps * s;
        const ChFrame<> frame = CarrierFrame(c);
        auto distance = [&](const ChVector3d& world) -> float {
            const ChVector3d p = frame.TransformPointParentToLocal(world);
            // Top of the load above the cavity's floor at this point
            const double r = std::hypot(p.x() - center.x(), p.y() - center.y());
            double top = level + std::max(0.0, cone - r * slope);
            top += lumps * (0.6 * Noise(p.x() / (3 * s), p.y() / (3 * s)) + 0.4 * Noise(p.x() / s + 17, p.y() / s + 5));
            const double above = p.z() - (c.cavity.min.z() + top);
            // Inside the cavity's walls
            const double out_x = std::abs(p.x() - center.x()) - 0.5 * size.x();
            const double out_y = std::abs(p.y() - center.y()) - 0.5 * size.y();
            const double below = c.cavity.min.z() - p.z();
            return static_cast<float>(std::max({above, out_x, out_y, below}));
        };
        const double reach = cone + lumps + 2 * s;
        ChAABB local(c.cavity.min - ChVector3d(s), ChVector3d(c.cavity.max.x() + s, c.cavity.max.y() + s, c.cavity.min.z() + level + reach));
        const ChAABB box = local.Transform(frame);
        m_grid.SplatFunction(box.min, box.max, distance);
    }

    // Particles out of carriers, as clods
    for (const auto& p : m_free)
        m_grid.SplatSphere(p.pos, p.size * s, blend);
    m_grid.Remesh();

    auto mesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    m_grid.AppendMeshes(*mesh);
    m_mesh = mesh;
    ++m_mesh_version;
    if (m_use_visual_model) {
        auto model = m_body->GetVisualModel();
        if (!model) {
            model = chrono_types::make_shared<ChVisualModel>();
            m_body->AddVisualModel(model);
        }
        model->Clear();
        if (mesh->GetNumTriangles() > 0) {
            auto shape = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
            shape->SetMesh(mesh);
            shape->SetMutable(false);
            shape->AddMaterial(m_material);
            model->AddShape(shape);
        }
    }
}

}  // namespace planet
}  // namespace chrono
