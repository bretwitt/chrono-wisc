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
//
// SCM soil stiffness that varies with the state of the wheel pressing on it,
// from a sparse model fitted to a reference soil model.
//
// =============================================================================

#ifndef SCM_SINDY_STIFFNESS_H
#define SCM_SINDY_STIFFNESS_H

#include <memory>
#include <unordered_map>

#include "chrono_vehicle/ChApiVehicle.h"
#include "chrono_vehicle/terrain/SCMSindyResidual.h"
#include "chrono_vehicle/terrain/SCMTerrain.h"

namespace chrono {
namespace vehicle {

/// @addtogroup vehicle_terrain
/// @{

/// SCM soil parameters with the Bekker Kphi under each wheel scaled by exp(f(state)), f a sparse model of the wheel's
/// state (the features of SCMSindyResidual: slip, slope, sinkage, rut ahead, slip history).
///
/// A residual on the normal force cannot match a stiff reference soil: at one sinkage it carries whatever load the
/// wheel has, so no function of the wheel's state gives its normal force. Scaling SCM's stiffness instead leaves the
/// load balance to SCM's own pressure-sinkage law, and moves where it settles: a wheel sinks as far into SCM as into
/// the reference soil. The map f is fitted to the multipliers that make SCM, driven through a reference-model test,
/// carry the reference model's normal force.
///
/// The map is an SCMSindyResidual whose normal-force output (index 2) is read as ln(Kphi multiplier); its other
/// outputs are not used, and it is not registered as a force correction. Each contact point takes the multiplier of
/// the wheel nearest it, within the wheel's radius plus a margin; points away from the wheels keep the base soil.
class CH_VEHICLE_API SCMSindyStiffness : public SCMTerrain::SoilParametersCallback {
  public:
    /// Base soil parameters, as SCMTerrain::SetSoilParameters takes them.
    struct Soil {
        double kphi, kc, n, cohesion, friction, janosi, elastic_k, damping;
    };

    /// Stiffness from `map` on the wheels of `terrain` over the `base` soil. Add the wheels to the map (AddWheel) and set
    /// its terrain; register this with the terrain (SCMTerrain::RegisterSoilParametersCallback).
    SCMSindyStiffness(std::shared_ptr<SCMSindyResidual> map, const SCMTerrain* terrain, const Soil& base);

    /// Range the multiplier is held to (default: 0.2 to 20).
    void SetMultiplierRange(double lo, double hi) {
        m_lo = lo;
        m_hi = hi;
    }

    /// Kphi multiplier on a wheel at the current time (1 if the wheel is off the ground).
    double GetMultiplier(const ChBody& wheel) const;

    /// SCMTerrain::SoilParametersCallback
    virtual void Set(const ChVector3d& loc,
                     double& Bekker_Kphi,
                     double& Bekker_Kc,
                     double& Bekker_n,
                     double& Mohr_cohesion,
                     double& Mohr_friction,
                     double& Janosi_shear,
                     double& elastic_K,
                     double& damping_R) override;

  private:
    struct Cached {
        double time = -1;
        double multiplier = 1;
        bool on_ground = false;
        double radius = 0;
    };
    const Cached& Update(const ChBody& wheel) const;

    std::shared_ptr<SCMSindyResidual> m_map;
    const SCMTerrain* m_terrain;
    Soil m_base;
    double m_lo = 0.2, m_hi = 20;
    mutable std::unordered_map<const ChBody*, Cached> m_cache;
};

/// @} vehicle_terrain

}  // namespace vehicle
}  // namespace chrono

#endif
