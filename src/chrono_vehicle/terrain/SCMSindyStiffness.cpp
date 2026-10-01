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

#include "chrono_vehicle/terrain/SCMSindyStiffness.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "chrono/physics/ChSystem.h"

namespace chrono {
namespace vehicle {

SCMSindyStiffness::SCMSindyStiffness(std::shared_ptr<SCMSindyResidual> map, const SCMTerrain* terrain, const Soil& base)
    : m_map(std::move(map)), m_terrain(terrain), m_base(base) {
    if (!m_map || !m_terrain)
        throw std::invalid_argument("SCMSindyStiffness: null map or terrain");
}

const SCMSindyStiffness::Cached& SCMSindyStiffness::Update(const ChBody& wheel) const {
    // Once per wheel per time: its state and multiplier
    const double now = wheel.GetSystem() ? wheel.GetSystem()->GetChTime() : 0.0;
    Cached& c = m_cache[&wheel];
    if (c.time == now)
        return c;
    c.time = now;
    c.multiplier = 1;
    SCMSindyResidual::State st;
    c.on_ground = m_map->GetState(wheel, VNULL, st);
    if (c.on_ground) {
        c.radius = st.radius;
        if (m_map->IsTrained())
            c.multiplier = std::clamp(std::exp(m_map->Evaluate(st.x)[2]), m_lo, m_hi);
    }
    return c;
}

double SCMSindyStiffness::GetMultiplier(const ChBody& wheel) const {
    return Update(wheel).multiplier;
}

void SCMSindyStiffness::Set(const ChVector3d& loc,
                            double& Bekker_Kphi,
                            double& Bekker_Kc,
                            double& Bekker_n,
                            double& Mohr_cohesion,
                            double& Mohr_friction,
                            double& Janosi_shear,
                            double& elastic_K,
                            double& damping_R) {
    Bekker_Kc = m_base.kc;
    Bekker_n = m_base.n;
    Mohr_cohesion = m_base.cohesion;
    Mohr_friction = m_base.friction;
    Janosi_shear = m_base.janosi;
    damping_R = m_base.damping;

    // The wheel nearest the contact point, in SCM's plane
    const ChFrame<> frame(m_terrain->GetReferenceFrame());
    double multiplier = 1;
    double best = 1e300;
    for (ChBody* wheel : m_map->GetBodies()) {
        const ChVector3d p = frame.TransformPointParentToLocal(wheel->GetFrameRefToAbs().GetPos());
        const double d2 = (p.x() - loc.x()) * (p.x() - loc.x()) + (p.y() - loc.y()) * (p.y() - loc.y());
        const Cached& c = Update(*wheel);
        if (c.on_ground && d2 < best && d2 < std::pow(c.radius + 0.1, 2)) {
            best = d2;
            multiplier = c.multiplier;
        }
    }
    Bekker_Kphi = m_base.kphi * multiplier;
    elastic_K = std::max(m_base.elastic_k, 5 * Bekker_Kphi);
}

}  // namespace vehicle
}  // namespace chrono
