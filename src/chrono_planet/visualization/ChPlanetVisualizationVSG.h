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
// VSG plugin drawing the quadtree terrain around the camera in a site frame.
//
// =============================================================================

#ifndef CH_PLANET_VISUALIZATION_VSG_H
#define CH_PLANET_VISUALIZATION_VSG_H

#include <cstdint>
#include <memory>
#include <unordered_map>

#include "chrono/assets/ChVisualMaterial.h"

#include "chrono_vsg/ChVisualSystemVSG.h"

#include "chrono_planet/ChApiPlanet.h"
#include "chrono_planet/ChSiteFrame.h"
#include "chrono_planet/lod/QuadtreeWorld.h"

namespace chrono {
namespace planet {

/// @addtogroup planet_module
/// @{

/// VSG plugin that mirrors the quadtree's resident tiles into the scene.
/// Each frame the level of detail is stepped at the camera position on simulation time, and every
/// tile mesh the world holds is re-projected through the site frame, so the drawn ground is the
/// surface PlanetTerrain and PlanetSCMTerrain put under the wheels. Tiles use the visual system's
/// PBR pipeline with one material; this is the geometry path, ahead of the Hapke shading port.
class CH_PLANET_API ChPlanetVisualizationVSG : public vsg3d::ChVisualSystemVSGPlugin {
  public:
    /// Construct the plugin over a quadtree world, drawn in the given site frame.
    ChPlanetVisualizationVSG(std::shared_ptr<QuadtreeWorld> world, const ChSiteFrame& site);
    ~ChPlanetVisualizationVSG();

    /// Set the material every tile is drawn with (default: a matte regolith grey).
    void SetMaterial(std::shared_ptr<ChVisualMaterial> material) { m_material = material; }

    /// Draw tiles as wireframe (default: false). Applies to tiles built after the call.
    void SetWireframe(bool val) { m_wireframe = val; }

    /// Show or hide the terrain.
    void SetVisible(bool val);

    /// Number of tiles currently in the scene.
    size_t GetNumTiles() const { return m_tiles.size(); }

    /// The quadtree world this plugin drives.
    std::shared_ptr<QuadtreeWorld> GetWorld() const { return m_world; }

    virtual void OnAttach() override;
    virtual void OnBindAssets() override;
    virtual void OnRender() override;

  private:
    using NodeKey = const QuadTree<TileMetadata, Spherical>*;

    struct Tile {
        std::uint64_t mesh_id;
        vsg::ref_ptr<vsg::Node> node;
    };

    /// Build a VSG subgraph for one tile mesh, in site coordinates.
    vsg::ref_ptr<vsg::Node> BuildTile(const Mesh& mesh) const;

    /// Add and drop tiles so the scene matches the world's resident meshes.
    void SyncTiles();

    std::shared_ptr<QuadtreeWorld> m_world;
    ChSiteFrame m_site;
    std::shared_ptr<ChVisualMaterial> m_material;
    bool m_wireframe;
    bool m_visible;

    vsg::ref_ptr<vsg::Switch> m_terrain_scene;  ///< VSG scene holding the tiles
    std::unordered_map<NodeKey, Tile> m_tiles;  ///< tiles in the scene, by quadtree node
    unsigned long long m_mesh_version;          ///< world mesh set version last mirrored
};

/// @} planet_module

}  // namespace planet
}  // namespace chrono

#endif
