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
// A body that digs a ChSiteVolume: it cuts the soil it sweeps through, carries
// it, dumps it, and feels the soil's resistance.
//
// =============================================================================

#ifndef CH_EXCAVATION_TOOL_H
#define CH_EXCAVATION_TOOL_H

#include <memory>

#include "chrono/core/ChFrame.h"
#include "chrono/physics/ChBody.h"
#include "chrono/utils/ChConstants.h"

#include "chrono_planet/ChApiPlanet.h"
#include "chrono_planet/volume/ChSdfShape.h"
#include "chrono_planet/volume/ChSiteVolume.h"
#include "chrono_planet/volume/ChSoilParticles.h"

namespace chrono {
namespace planet {

/// @addtogroup planet_module
/// @{

/// An excavation tool, such as a bucket, a blade or a drum's scoop, on a body: at each Update, the soil of a
/// ChSiteVolume that its cutting shape swept through since the last edit is removed and carried as payload, up
/// to a capacity; soil cut beyond it is pushed ahead of the tool as a pile. Dump pours the payload onto the
/// ground. Soil volumes are in place (bank) volumes, so what is dug and dumped balances; a swell factor makes
/// dumped soil take more room.
///
/// With particles (SetParticles), the soil the tool carries stays in view: as it is cut it fills the tool's cavity
/// as particles, and when the tool's open top tips past the angle of repose it pours out, faster the more it tips,
/// as particles that fall and land back in the ground. Dump then lets the whole load fall from where the tool is.
///
/// The soil resists the tool with the fundamental equation of earthmoving (Reece 1964; McKyes 1985, Soil
/// Cutting and Tillage): a blade of width w and rake angle beta cutting at depth d needs a force
///   P = w (gamma g d^2 N_gamma + c d N_c + q d N_q),
/// the N factors given by the soil's friction angle phi, the soil-tool friction angle delta and the angle rho of
/// the failure plane in front of the blade, which is chosen to make P least. P acts opposite to the tool's
/// horizontal motion and presses it down, at the angle beta + delta from the vertical. The payload's weight is
/// added. The depth is how far the lowest point of the cutting shape is below the top of the ground, and the
/// width the shape's extent across the motion. The forces go into an accumulator of the body.
class CH_PLANET_API ChExcavationTool {
  public:
    /// Soil properties.
    struct Soil {
        double bulk_density = 1600.0;                    ///< in-place density (kg/m^3)
        double cohesion = 1000.0;                        ///< c (Pa)
        double friction_angle = 35 * CH_DEG_TO_RAD;       ///< phi (rad)
        double tool_friction_angle = 20 * CH_DEG_TO_RAD;  ///< delta, between soil and tool (rad)
        double surcharge = 0.0;                          ///< q, pressure on the ground in front of the blade (Pa)
        double repose_angle = 35 * CH_DEG_TO_RAD;         ///< slope of piles (rad)
        double swell = 1.0;                              ///< loose volume of dumped soil over its bank volume
        double gravity = 1.62;                           ///< gravitational acceleration (m/s^2), along -z
    };

    /// A tool cutting with `cutter`, placed on the body by `cutter_frame` (in the body's reference frame).
    ChExcavationTool(std::shared_ptr<ChSiteVolume> volume,
                     std::shared_ptr<ChBody> body,
                     std::shared_ptr<ChSdfShape> cutter,
                     const ChFrame<>& cutter_frame = ChFrame<>());

    void SetSoil(const Soil& soil) { m_soil = soil; }
    const Soil& GetSoil() const { return m_soil; }

    /// Bank volume of soil the tool can carry (m^3, default: unlimited).
    void SetCapacity(double volume) { m_capacity = volume; }
    double GetCapacity() const { return m_capacity; }

    /// Angle of the cutting face to the horizontal (rad, default: 45 degrees).
    void SetRakeAngle(double angle) { m_rake = angle; }
    double GetRakeAngle() const { return m_rake; }

    /// Show the payload as particles held in `cavity`, a box in the cutting shape's frame whose open side is +z, and
    /// pour it by tipping the tool (see the class description). Call before the first Update.
    void SetParticles(std::shared_ptr<ChSoilParticles> particles, const ChAABB& cavity);
    std::shared_ptr<ChSoilParticles> GetParticles() const { return m_particles; }

    /// Apply the soil's resistance and the payload's weight to the body (default: true).
    void EnableForces(bool val) { m_forces = val; }

    /// Cut the soil the tool swept since the last edit, and apply the soil's forces. Call once per step, before
    /// advancing the system. The first call only records where the tool is.
    void Update();

    /// Pour the payload onto the ground below the cutting shape, or with particles, let it fall from the tool. Returns
    /// the bank volume poured (m^3); with particles, less than a particle's worth stays in the tool.
    double Dump();

    /// Pour the payload onto the ground at site x/y at once, without particles. Returns the bank volume poured (m^3).
    double Dump(double x, double y);

    double GetPayloadVolume() const { return m_payload; }                            ///< carried (m^3)
    double GetPayloadMass() const { return m_payload * m_soil.bulk_density; }        ///< carried (kg)
    /// Cut in all (m^3), counting soil pushed ahead of a full tool again each time the tool cuts into it.
    double GetExcavatedVolume() const { return m_excavated; }
    double GetDumpedVolume() const { return m_dumped; }                              ///< dumped or poured (m^3)
    double GetSpilledVolume() const { return m_spilled; }                            ///< pushed ahead (m^3)
    /// Cut beyond the capacity and not yet piled ahead (m^3): piled once enough has gathered.
    double GetPendingSpillVolume() const { return m_spill_pending; }
    double GetCuttingDepth() const { return m_depth; }                               ///< at the last Update (m)
    const ChVector3d& GetResistanceForce() const { return m_force; }                 ///< at the last Update (N)

    /// Force (N) the fundamental equation of earthmoving gives for a blade of width w (m) and rake angle beta (rad)
    /// at depth d (m) in the given soil, for the failure plane angle that makes it least.
    static double CuttingForce(const Soil& soil, double width, double depth, double beta);

    std::shared_ptr<ChBody> GetBody() const { return m_body; }
    std::shared_ptr<ChSiteVolume> GetVolume() const { return m_volume; }

  private:
    ChFrame<> CutterFrame() const;
    double Pour(double x, double y, double volume);

    std::shared_ptr<ChSiteVolume> m_volume;
    std::shared_ptr<ChBody> m_body;
    std::shared_ptr<ChSdfShape> m_cutter;
    ChFrame<> m_cutter_frame;
    Soil m_soil;
    double m_capacity = 1e300;
    double m_rake = 45 * CH_DEG_TO_RAD;
    bool m_forces = true;
    unsigned int m_accumulator;

    bool m_started = false;
    ChFrame<> m_last;  // cutter placement at the last edit
    double m_last_time = 0;

    std::shared_ptr<ChSoilParticles> m_particles;
    int m_carrier = -1;
    double m_pour_load = 0;     // payload when the tool tipped past the angle of repose, 0 when upright
    double m_pour_pending = 0;  // volume due to pour, less than a particle

    double m_payload = 0;
    double m_excavated = 0;
    double m_dumped = 0;
    double m_spilled = 0;
    double m_spill_pending = 0;
    double m_depth = 0;
    ChVector3d m_force;
};

/// @} planet_module

}  // namespace planet
}  // namespace chrono

#endif
