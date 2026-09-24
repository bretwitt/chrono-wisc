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
// VSG plugin drawing the soil particles of a ChSoilParticles: carried, poured,
// falling and landed.
//
// =============================================================================

#ifndef CH_SOIL_PARTICLES_VISUALIZATION_VSG_H
#define CH_SOIL_PARTICLES_VISUALIZATION_VSG_H

#include <cstdint>
#include <memory>

#include "chrono_vsg/ChVisualSystemVSG.h"

#include "chrono_planet/ChApiPlanet.h"
#include "chrono_planet/visualization/ChMeshDrawsVSG.h"
#include "chrono_planet/volume/ChSoilParticles.h"

namespace chrono {
namespace planet {

/// @addtogroup planet_module
/// @{

/// VSG plugin that draws the mesh of a ChSoilParticles as it changes. Turn the particles' visual model off
/// (ChSoilParticles::SetUseVisualModel) and call ChSoilParticles::UpdateMesh before rendering.
class CH_PLANET_API ChSoilParticlesVisualizationVSG : public vsg3d::ChVisualSystemVSGPlugin {
  public:
    explicit ChSoilParticlesVisualizationVSG(std::shared_ptr<ChSoilParticles> particles);
    ~ChSoilParticlesVisualizationVSG();

    /// Show or hide the particles.
    void SetVisible(bool val);

    virtual void OnBindAssets() override;
    virtual void OnRender() override;

  private:
    std::shared_ptr<ChSoilParticles> m_particles;
    vsg::ref_ptr<vsg::Switch> m_scene;
    std::unique_ptr<ChMeshDrawsVSG> m_draws;
    ChMeshDrawsVSG::Draws m_current;
    std::uint64_t m_version = ~0ull;
    bool m_visible = true;
};

/// @} planet_module

}  // namespace planet
}  // namespace chrono

#endif
