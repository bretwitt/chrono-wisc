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
// Soil out of the ground of a ChSiteVolume and on its way back: carried in a
// tool, poured, falling, and landing, as particles drawn as one smooth mass.
//
// =============================================================================

#ifndef CH_SOIL_PARTICLES_H
#define CH_SOIL_PARTICLES_H

#include <cstdint>
#include <memory>
#include <random>
#include <vector>

#include "chrono/assets/ChVisualMaterial.h"
#include "chrono/assets/ChVisualShapeTriangleMesh.h"
#include "chrono/core/ChFrame.h"
#include "chrono/geometry/ChAABB.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChSystem.h"
#include "chrono/utils/ChConstants.h"

#include "chrono_planet/ChApiPlanet.h"
#include "chrono_planet/volume/ChSdfMesher.h"
#include "chrono_planet/volume/ChSiteVolume.h"

namespace chrono {
namespace planet {

/// @addtogroup planet_module
/// @{

/// Soil between leaving the ground of a ChSiteVolume and going back into it, as particles, so it stays in view the
/// whole way: carried in a tool, pouring out of it, falling, and landing. Each particle stands for a cube of bank
/// soil `spacing` on a side. A carrier's load is drawn as the soil filling its cavity: level at the height that holds
/// the load, heaped at the angle of repose above the rim, its surface lumpy with clods. Particles out of a carrier are
/// drawn as clods, spheres of a little varied size barely blended, so a pour reads as falling grains. All of it is
/// splatted into one sparse lattice of bricks (ChSparseSdfGrid), which exist only where soil is and are meshed again
/// only when what is in them changes: soil at rest costs nothing to draw again.
///
/// Carriers, such as the bucket of a ChExcavationTool (see ChExcavationTool::SetParticles), hold particles in a
/// lattice over a cavity in their own frame, filled from the bottom up. Release lets particles go, the lowest first,
/// at the carrier's velocity; they fall under gravity and rest where they meet the ground. Resting particles are
/// merged into the ground every `merge_period`, those within a cell `merge_cell` wide poured as one pile of their
/// volume at the angle of repose (ChSiteVolume::Deposit), and are drawn until then.
///
/// Call Advance once per step and UpdateMesh before drawing. The mesh is in site coordinates; draw it with
/// ChSoilParticlesVisualizationVSG, or through the visual model of GetBody() for renderers that read the system's
/// visual models, such as Chrono::Sensor. Not thread-safe.
class CH_PLANET_API ChSoilParticles {
  public:
    struct Params {
        double spacing = 0.04;                     ///< edge of the cube of bank soil a particle stands for (m)
        double gravity = 1.62;                     ///< gravitational acceleration (m/s^2), along -z
        double repose_angle = 35 * CH_DEG_TO_RAD;  ///< slope of the piles landed particles make (rad)
        double pour_time = 1.0;                    ///< time a load takes to pour from a tool tipped past vertical (s)
        double merge_period = 0.05;                ///< landed particles merge into the ground this often (s)
        double merge_cell = 0.15;                  ///< landed particles within a cell this wide merge as one pile (m)
        double blend = 0.15;                       ///< width of the smooth union of clods, in spacings
        double lumps = 0.25;                       ///< height of the lumps on a carried load, in spacings
        double spread = 0.12;                      ///< random speed particles leave a carrier with, besides its own (m/s)
    };

    /// Particles going to and from the ground of `volume`. A fixed body for their visual shape is added to `sys`.
    ChSoilParticles(ChSystem* sys, std::shared_ptr<ChSiteVolume> volume, const Params& params);
    /// Particles with the default parameters.
    ChSoilParticles(ChSystem* sys, std::shared_ptr<ChSiteVolume> volume);
    ~ChSoilParticles();

    const Params& GetParams() const { return m_params; }
    /// Bank volume one particle stands for (m^3).
    double GetParticleVolume() const { return m_params.spacing * m_params.spacing * m_params.spacing; }

    /// Material the particles are drawn with (default: a darker regolith gray).
    void SetMaterial(std::shared_ptr<ChVisualMaterial> material) { m_material = material; }
    std::shared_ptr<ChVisualMaterial> GetMaterial() const { return m_material; }

    /// Keep the particles' mesh as a visual shape of GetBody() (default: true), for renderers that read the system's
    /// visual models. Turn it off when ChSoilParticlesVisualizationVSG draws them.
    void SetUseVisualModel(bool val) { m_use_visual_model = val; }

    // -------------------------------------------------------------------------
    // Carriers

    /// Add a carrier: a cavity (a box in `frame`, which is in the body's reference frame) that holds particles.
    /// Returns its id.
    int AddCarrier(std::shared_ptr<ChBody> body, const ChFrame<>& frame, const ChAABB& cavity);

    /// Make a carrier hold `volume` of bank soil (m^3), in whole particles, adding or dropping the latest added.
    void SetCarried(int carrier, double volume);

    /// Let up to `volume` (m^3) of a carrier's particles go, the lowest first, at the carrier's velocity where they
    /// are plus a random `spread`, and from a little off their places in the load, so a pour breaks into a stream of
    /// clods. Returns the volume released, in whole particles.
    double Release(int carrier, double volume);

    /// Bank volume a carrier holds (m^3).
    double GetCarriedVolume(int carrier) const;

    // -------------------------------------------------------------------------
    // Motion and drawing

    /// Move falling particles on by dt (s), rest those that meet the ground, and merge resting ones into it when due.
    void Advance(double dt);

    /// Merge every resting particle into the ground now.
    void Merge();

    /// Mesh the particles, carried and free, in site coordinates.
    void UpdateMesh();

    /// The latest mesh (UpdateMesh), and a count that changes with it.
    std::shared_ptr<ChTriangleMeshConnected> GetMesh() const { return m_mesh; }
    std::uint64_t GetMeshVersion() const { return m_mesh_version; }

    /// The fixed body holding the particles' visual shape.
    std::shared_ptr<ChBody> GetBody() const { return m_body; }

    size_t GetNumCarried() const;
    size_t GetNumFalling() const;
    size_t GetNumResting() const;
    /// Bank volume of the free particles, falling or resting (m^3).
    double GetFreeVolume() const { return (GetNumFalling() + GetNumResting()) * GetParticleVolume(); }
    /// Bank volume merged into the ground (m^3), as the volume measured it.
    double GetMergedVolume() const { return m_merged; }
    /// Bank volume that did not go back into the ground (m^3): particles that landed beyond its edit region.
    double GetLostVolume() const { return m_lost; }

  private:
    struct Carrier {
        std::shared_ptr<ChBody> body;
        ChFrame<> frame;
        ChAABB cavity;
        int nx, ny;                 // lattice slots across the cavity; layers go up without bound
        std::vector<int> slots;     // slot of each carried particle, in the order added
    };
    struct Particle {
        ChVector3d pos, vel;
        bool resting;
        float size;  // radius of its clod, in spacings
    };

    ChFrame<> CarrierFrame(const Carrier& c) const { return c.body->GetFrameRefToAbs() * c.frame; }
    ChVector3d SlotPosition(const Carrier& c, int slot) const;

    ChSystem* m_system;
    std::shared_ptr<ChSiteVolume> m_volume;
    Params m_params;
    std::shared_ptr<ChBody> m_body;
    std::shared_ptr<ChVisualMaterial> m_material;
    bool m_use_visual_model = true;

    std::vector<Carrier> m_carriers;
    std::vector<Particle> m_free;
    std::mt19937 m_rng{7};
    double m_since_merge = 0;
    double m_merged = 0;
    double m_lost = 0;

    ChSparseSdfGrid m_grid;
    std::shared_ptr<ChTriangleMeshConnected> m_mesh;
    std::uint64_t m_mesh_version = 0;
};

/// @} planet_module

}  // namespace planet
}  // namespace chrono

#endif
