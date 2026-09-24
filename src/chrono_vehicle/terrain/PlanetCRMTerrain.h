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
// CRM soil over a window of a planet work site, taken from and given back to
// the site's signed distance field (planet::ChSiteVolume).
//
// =============================================================================

#ifndef PLANET_CRM_TERRAIN_H
#define PLANET_CRM_TERRAIN_H

#include <functional>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

#include "chrono_vehicle/ChApiVehicle.h"
#include "chrono_vehicle/terrain/CRMTerrain.h"
#include "chrono_fsi/sph/ChFsiFluidSystemSPH.h"

#include "chrono_planet/ChSiteFrame.h"
#include "chrono_planet/filters/ChDeformationFilter.h"
#include "chrono_planet/volume/ChSiteVolume.h"

namespace chrono {
namespace vehicle {

/// @addtogroup vehicle_terrain
/// @{

/// CRM terrain (SPH continuum soil, Chrono::FSI) over a rectangle of a planet site, taken from and given back to
/// ground the rest of the scene draws and collides with: a planet::ChSiteVolume at a work site (ConstructFromVolume,
/// PublishToVolume), or a height field with ruts in a planet::ChDeformationFilter where a rover drives
/// (ConstructFromHeight, PublishToDeformation; see PlanetCRMWindow for a window that follows the rover).
///
/// With a volume: The soil is taken from the volume: a particle at every lattice point, `spacing` apart,
/// where the volume holds soil, down to a floor below the lowest ground in the window, so trenches, undercuts and
/// piles dug before carry over. PublishToVolume gives the soil back: the particles are splatted into the volume over
/// the window, as the smooth union of spheres, so what the soil model moved shows in the volume's meshes and
/// collision shapes (planet::ChSiteVolumeShapes), and the change of soil volume is measured as for any edit.
///
/// Positions are in the site frame, z up. Set the soil and SPH parameters, add the bodies that touch the soil
/// (AddRigidBody), construct from the volume, and Initialize, as for any CRMTerrain. The window should lie within the
/// volume's edit region; around it the volume's own ground stays as it was.
class CH_VEHICLE_API PlanetCRMTerrain : public CRMTerrain {
  public:
    /// CRM soil of the given particle spacing (m), over ground held by `volume` if one is given.
    PlanetCRMTerrain(ChSystem& sys, double spacing, std::shared_ptr<planet::ChSiteVolume> volume = nullptr);
    ~PlanetCRMTerrain();

    /// Room above the highest ground in the window, and around the window, that the soil may be carried to (m,
    /// defaults: 2 and 1): the computational domain of the soil model. Particles that leave it drop out of the
    /// simulation, so it must hold what tools lift. Call before ConstructFromVolume.
    void SetHeadroom(double above, double around) {
        m_headroom = above;
        m_margin = around;
    }

    /// Fill `window` with particles where the volume holds soil, down to `depth` (m) below the lowest ground in it,
    /// with boundary markers under that floor and, for `side_flags` (fsi::sph::BoxSide), walls on the window's sides.
    void ConstructFromVolume(const planet::ChSiteRegion& window,
                             double depth,
                             int side_flags = fsi::sph::BoxSide::ALL & ~fsi::sph::BoxSide::Z_POS);

    /// Fill `window` with particles below the ground a height function gives (site x/y to z, m), down to `depth` (m)
    /// below the lowest ground in it, as ConstructFromVolume. Each column of particles is shifted up or down by less
    /// than half a spacing so the soil's surface lies at the ground, not below it by up to a spacing.
    void ConstructFromHeight(const std::function<double(double x, double y)>& height,
                             const planet::ChSiteRegion& window,
                             double depth,
                             int side_flags = fsi::sph::BoxSide::Z_NEG);

    /// State of a soil particle, to carry soil over from one soil model to another.
    struct ParticleState {
        ChVector3d pos;
        ChVector3d vel;
        double rho;
        double pressure;
        double mu;
        ChVector3d tau_diag;
        ChVector3d tau_offdiag;
    };

    /// The state of the soil's particles.
    std::vector<ParticleState> GetParticleStates() const;

    /// Soil to take over, as it is, where it lies over `region` (site x/y): those particles are seeded with their
    /// state in place of the lattice's, which then fills only the columns under them. Soil carried over from a window
    /// being replaced keeps the stress and compaction a rover has given it, so the rover does not sink into fresh soil
    /// when the window moves. Call before ConstructFromHeight.
    void SetCarriedParticles(std::vector<ParticleState> particles, const planet::ChSiteRegion& region);

    /// Seed the soil carrying its own weight: each particle starts at the pressure of the soil above it, up to the
    /// ground a height function gives (ConstructFromHeight), instead of at rest with no stress, which leaves the top
    /// layers unconfined and weak until they settle. Call before ConstructFromHeight (default: false).
    void SetInitialOverburden(bool val) { m_overburden = val; }

    /// Write how far the soil's top has moved since the window was seeded into a deformation filter, over the window,
    /// added to the height changes the filter held there then: ruts from before carry over, and soil the particles did
    /// not move shows no change. The first call, which should come right after Initialize, only records where the
    /// soil starts. Returns the nodes changed by more than `tolerance` (m).
    size_t PublishToDeformation(planet::ChDeformationFilter& filter, double tolerance = 5e-4);

    /// Write the particles into the volume over the window, above the floor: the ground there becomes the smooth
    /// union of spheres around them. Returns the change of soil volume in the volume (m^3).
    double PublishToVolume();

    const planet::ChSiteRegion& GetWindow() const { return m_window; }
    /// Height of the floor of the soil (m).
    double GetFloor() const { return m_floor; }
    std::shared_ptr<planet::ChSiteVolume> GetVolume() const { return m_volume; }

  private:
    // Initial state of the particles: carried over, or at rest (with the soil's weight above, if m_overburden)
    class SeedCallback : public fsi::sph::ChFsiFluidSystemSPH::ParticlePropertiesCallback {
      public:
        explicit SeedCallback(const PlanetCRMTerrain& terrain) : m_terrain(terrain) {}
        virtual void set(const fsi::sph::ChFsiFluidSystemSPH& sysSPH, const ChVector3d& pos) override;

      private:
        const PlanetCRMTerrain& m_terrain;
    };

    virtual ChVector3d Grid2Point(const ChVector3i& p) override;
    static int64_t CarriedKey(const ChVector3d& pos);

    void Build(const planet::ChSiteRegion& window,
               double depth,
               int side_flags,
               const std::function<double(double, double)>& top_at,
               const std::function<bool(const ChVector3d&)>& soil_at);

    std::shared_ptr<planet::ChSiteVolume> m_volume;
    int m_nx = 0, m_ny = 0;           // lattice columns over the window
    int m_k_floor = 0;                // lattice level of the floor
    std::vector<double> m_shift;      // height shift of each column's particles (ConstructFromHeight)
    std::vector<double> m_column_ground;     // ground over each column (ConstructFromHeight)
    bool m_overburden = false;
    std::vector<ParticleState> m_carried;               // soil taken over as it is (SetCarriedParticles)
    planet::ChSiteRegion m_carried_region;
    std::unordered_map<int64_t, size_t> m_carried_index;  // carried particles by quantized position
    static constexpr int kCarried = -(1 << 29);          // lattice x index marking a carried particle (y: its index)
    std::vector<float> m_reference;   // top of the soil over each column when first published
    std::vector<float> m_base;        // height changes of the filter's nodes over the window then
    planet::ChSiteRegion m_window;
    double m_floor = 0;
    double m_headroom = 2.0;
    double m_margin = 1.0;
};

/// @} vehicle_terrain

}  // namespace vehicle
}  // namespace chrono

#endif
