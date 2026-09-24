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

#include "chrono_planet/volume/ChExcavationTool.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace chrono {
namespace planet {

namespace {
// Soil pushed ahead of a full tool is piled once this much has gathered (m^3)
constexpr double kSpillBatch = 1e-3;
// Speed (m/s) over which the resistance ramps in, so a tool at rest feels none
constexpr double kRestSpeed = 0.01;

double Cot(double a) {
    return 1.0 / std::tan(a);
}
}  // namespace

ChExcavationTool::ChExcavationTool(std::shared_ptr<ChSiteVolume> volume,
                                   std::shared_ptr<ChBody> body,
                                   std::shared_ptr<ChSdfShape> cutter,
                                   const ChFrame<>& cutter_frame)
    : m_volume(std::move(volume)), m_body(std::move(body)), m_cutter(std::move(cutter)), m_cutter_frame(cutter_frame), m_force(VNULL) {
    if (!m_volume || !m_body || !m_cutter)
        throw std::invalid_argument("ChExcavationTool: null volume, body or cutter");
    m_accumulator = m_body->AddAccumulator();
}

void ChExcavationTool::SetParticles(std::shared_ptr<ChSoilParticles> particles, const ChAABB& cavity) {
    m_particles = std::move(particles);
    m_carrier = m_particles ? m_particles->AddCarrier(m_body, m_cutter_frame, cavity) : -1;
}

ChFrame<> ChExcavationTool::CutterFrame() const {
    return m_body->GetFrameRefToAbs() * m_cutter_frame;
}

double ChExcavationTool::CuttingForce(const Soil& soil, double width, double depth, double beta) {
    if (!(width > 0) || !(depth > 0))
        return 0;
    const double phi = soil.friction_angle, delta = soil.tool_friction_angle;
    const double gamma_g = soil.bulk_density * soil.gravity;
    beta = std::clamp(beta, 5 * CH_DEG_TO_RAD, 175 * CH_DEG_TO_RAD);
    double best = 1e300;
    for (int n = 1; n < 90; ++n) {
        const double rho = n * CH_DEG_TO_RAD;
        const double denom = std::cos(beta + delta) + std::sin(beta + delta) * Cot(rho + phi);
        if (!(denom > 1e-6))
            continue;
        const double n_gamma = (Cot(beta) + Cot(rho)) / (2 * denom);
        const double n_c = (1 + Cot(rho) * Cot(rho + phi)) / denom;
        const double n_q = (Cot(beta) + Cot(rho)) / denom;
        const double force = width * (gamma_g * depth * depth * n_gamma + soil.cohesion * depth * n_c + soil.surcharge * depth * n_q);
        if (force > 0)
            best = std::min(best, force);
    }
    return best < 1e300 ? best : 0;
}

void ChExcavationTool::Update() {
    const ChFrame<> now = CutterFrame();
    // The system's clock: a body's own is only brought up to date when the system updates it
    const double time = m_body->GetSystem() ? m_body->GetSystem()->GetChTime() : m_body->GetChTime();
    m_body->EmptyAccumulator(m_accumulator);
    m_force = VNULL;
    m_depth = 0;
    if (!m_started) {
        m_last = now;
        m_last_time = time;
        m_started = true;
        return;
    }
    const double dt = std::max(0.0, time - m_last_time);
    m_last_time = time;

    // How the tool moves: point speeds are taken in the body's center of mass frame
    const ChVector3d center = now.GetPos();
    const ChVector3d velocity = m_body->PointSpeedLocalToParent(m_body->TransformPointParentToLocal(center));
    ChVector3d ahead(velocity.x(), velocity.y(), 0);
    const double speed = ahead.Length();
    if (speed > 1e-9)
        ahead /= speed;

    // The leading bottom of the cutting shape, and how deep it is in the ground ahead of it, before this cut
    const ChAABB box = m_cutter->GetBoundingBox();
    ChVector3d corners[8];
    double travel = 0, bottom = 1e300;
    for (int c = 0; c < 8; ++c) {
        const ChVector3d corner((c & 1) ? box.max.x() : box.min.x(), (c & 2) ? box.max.y() : box.min.y(), (c & 4) ? box.max.z() : box.min.z());
        corners[c] = now.TransformPointLocalToParent(corner);
        travel = std::max(travel, (corners[c] - m_last.TransformPointLocalToParent(corner)).Length());
        bottom = std::min(bottom, corners[c].z());
    }
    const double h = m_volume->GetVoxelSize();
    ChVector3d edge = center;
    double lead = -1e300;
    for (const auto& p : corners)
        if (p.z() < bottom + h && p.Dot(ahead) > lead) {
            lead = p.Dot(ahead);
            edge = p;
        }
    const ChVector3d probe = edge + ahead * h;
    m_depth = std::max(0.0, m_volume->GetTopHeight(probe.x(), probe.y()) - edge.z());

    // Cut once the tool has moved a quarter voxel, so slow motions are not lost to the lattice
    if (travel > 0.25 * h) {
        const double removed = m_volume->SubtractSwept(*m_cutter, m_last, now);
        m_last = now;
        m_excavated += removed;
        const double take = std::clamp(m_capacity - m_payload, 0.0, removed);
        m_payload += take;
        m_spill_pending += removed - take;
    }

    // The load as particles, and poured while the tool is tipped past the angle of repose
    if (m_particles) {
        const auto& params = m_particles->GetParams();
        const double tilt = std::acos(std::clamp(now.GetRot().GetAxisZ().z(), -1.0, 1.0));
        if (tilt > params.repose_angle && m_payload > 0) {
            if (m_pour_load <= 0)
                m_pour_load = m_payload;
            const double rate = std::clamp((tilt - params.repose_angle) / (CH_PI_2 - params.repose_angle), 0.0, 1.0);
            m_pour_pending += m_pour_load / std::max(1e-6, params.pour_time) * rate * dt;
            const double poured = m_particles->Release(m_carrier, std::min(m_pour_pending, m_payload));
            m_pour_pending = std::max(0.0, m_pour_pending - poured);
            m_payload -= poured;
            m_dumped += poured;
        } else if (tilt <= params.repose_angle) {
            m_pour_load = 0;
            m_pour_pending = 0;
        }
        m_particles->SetCarried(m_carrier, m_payload);
    }

    // Soil the full tool cannot hold is pushed ahead of it, piled clear of the cutting shape so the tool does not
    // cut it again at once: ahead of its motion, or of its forward axis while it barely moves across the ground
    if (m_spill_pending > kSpillBatch) {
        ChVector3d dir = ahead;
        if (speed < kRestSpeed) {
            dir = now.GetRot().GetAxisX();
            dir.z() = 0;
            dir = dir.Length() > 1e-9 ? dir.GetNormalized() : ChVector3d(1, 0, 0);
        }
        double extent = 0;
        for (const auto& p : corners)
            extent = std::max(extent, (p - center).Dot(dir));
        const double pile_reach = std::cbrt(3 * m_spill_pending * m_soil.swell / (CH_PI * std::tan(m_soil.repose_angle)));
        const ChVector3d at = center + dir * (extent + pile_reach + h);
        Pour(at.x(), at.y(), m_spill_pending);
        m_spilled += m_spill_pending;
        m_spill_pending = 0;
    }

    if (!m_forces)
        return;

    // The soil's resistance, for the depth and width the tool cuts at, ramped in with its speed
    if (m_depth > 0 && speed > 1e-9) {
        const ChVector3d across(-ahead.y(), ahead.x(), 0);
        double lo = 1e300, hi = -1e300;
        for (const auto& p : corners) {
            lo = std::min(lo, p.Dot(across));
            hi = std::max(hi, p.Dot(across));
        }
        const double width = hi - lo;
        const double push = CuttingForce(m_soil, width, m_depth, m_rake);
        const double ramp = std::min(1.0, speed / kRestSpeed);
        const double angle = m_rake + m_soil.tool_friction_angle;
        m_force = (-ahead * std::sin(angle) - ChVector3d(0, 0, 1) * std::cos(angle)) * (push * ramp);
    }
    m_force += ChVector3d(0, 0, -GetPayloadMass() * m_soil.gravity);
    m_body->AccumulateForce(m_accumulator, m_force, edge, false);
}

double ChExcavationTool::Pour(double x, double y, double volume) {
    if (!(volume > 0))
        return 0;
    return m_volume->Deposit(x, y, volume * m_soil.swell, m_soil.repose_angle);
}

double ChExcavationTool::Dump() {
    if (m_particles) {
        const double released = m_particles->Release(m_carrier, m_payload);
        m_payload -= released;
        m_dumped += released;
        m_particles->SetCarried(m_carrier, m_payload);
        return released;
    }
    const ChVector3d at = CutterFrame().GetPos();
    return Dump(at.x(), at.y());
}

double ChExcavationTool::Dump(double x, double y) {
    const double volume = m_payload;
    m_payload = 0;
    if (m_particles)
        m_particles->SetCarried(m_carrier, 0);
    Pour(x, y, volume);
    m_dumped += volume;
    return volume;
}

}  // namespace planet
}  // namespace chrono
