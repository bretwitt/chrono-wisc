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
// VSG plugin drawing the ground of a ChSiteVolume as it is dug.
//
// =============================================================================

#ifndef CH_SITE_VOLUME_VISUALIZATION_VSG_H
#define CH_SITE_VOLUME_VISUALIZATION_VSG_H

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

#include "chrono_vsg/ChVisualSystemVSG.h"

#include "chrono_planet/ChApiPlanet.h"
#include "chrono_planet/visualization/ChMeshDrawsVSG.h"
#include "chrono_planet/volume/ChSiteVolumeShapes.h"

namespace chrono {
namespace planet {

/// @addtogroup planet_module
/// @{

/// VSG plugin that draws the brick meshes of a ChSiteVolumeShapes, adding and dropping them as edits rebuild
/// bricks. The VSG visual system binds a body's visual shapes once, so it cannot follow ground that changes;
/// turn the shapes' visual model off (ChSiteVolumeShapes::SetUseVisualModel) and attach this plugin instead.
/// Call ChSiteVolumeShapes::Update before rendering.
///
/// All bricks of a material share one set of pipeline state, built once, so a rebuilt brick costs only the
/// upload of its geometry: digging rebuilds tens of bricks a frame.
class CH_PLANET_API ChSiteVolumeVisualizationVSG : public vsg3d::ChVisualSystemVSGPlugin {
  public:
    explicit ChSiteVolumeVisualizationVSG(std::shared_ptr<ChSiteVolumeShapes> shapes);
    ~ChSiteVolumeVisualizationVSG();

    /// Draw the ground as wireframe (default: false). The bricks in the scene are rebuilt on the next frame.
    void SetWireframe(bool val) { m_wireframe_requested = val; }

    /// Show or hide the ground.
    void SetVisible(bool val);

    /// Number of brick meshes in the scene.
    size_t GetNumBricks() const { return m_bricks.size(); }

    virtual void OnBindAssets() override;
    virtual void OnRender() override;

  private:
    void Sync();

    struct BrickNodes {
        std::shared_ptr<ChVisualShapeTriangleMesh> shape;  ///< held, so its address is not reused
        ChMeshDrawsVSG::Draws draws;
    };

    std::shared_ptr<ChSiteVolumeShapes> m_shapes;
    vsg::ref_ptr<vsg::Switch> m_scene;
    std::unique_ptr<ChMeshDrawsVSG> m_draws;
    std::unordered_map<const ChVisualShapeTriangleMesh*, BrickNodes> m_bricks;
    std::uint64_t m_version = ~0ull;
    bool m_visible = true;
    bool m_wireframe = false;
    bool m_wireframe_requested = false;
};

/// @} planet_module

}  // namespace planet
}  // namespace chrono

#endif
