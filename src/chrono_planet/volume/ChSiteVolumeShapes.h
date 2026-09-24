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
// The ground of a ChSiteVolume as Chrono visual shapes and collision meshes,
// rebuilt brick by brick as it is dug.
//
// =============================================================================

#ifndef CH_SITE_VOLUME_SHAPES_H
#define CH_SITE_VOLUME_SHAPES_H

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

#include "chrono/assets/ChVisualMaterial.h"
#include "chrono/assets/ChVisualShapeTriangleMesh.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChContactMaterial.h"
#include "chrono/physics/ChSystem.h"

#include "chrono_planet/ChApiPlanet.h"
#include "chrono_planet/volume/ChSiteVolume.h"

namespace chrono {
namespace planet {

/// @addtogroup planet_module
/// @{

/// Draws the ground of a ChSiteVolume and, optionally, lets bodies collide with it. Each brick of the volume
/// holding ground gets a triangle mesh (ChSiteVolume::MeshBrick), as a visual shape on one fixed body and as a
/// collision mesh on a fixed body of its own, so an edit rebuilds only the bricks it changed. Call Update after
/// editing, for example once per rendered frame.
///
/// The brick bodies share a collision family that does not collide with itself. Other fixed bodies are still
/// tested against them, so give fixed bodies near the volume, such as a terrain patch, a family that does not
/// collide with GetCollisionFamily().
class CH_PLANET_API ChSiteVolumeShapes {
  public:
    /// Add the ground's bodies to `sys`.
    ChSiteVolumeShapes(ChSystem* sys, std::shared_ptr<ChSiteVolume> volume);
    ~ChSiteVolumeShapes();

    /// Material of the ground (default: gray). Applies to bricks meshed after the call.
    void SetMaterial(std::shared_ptr<ChVisualMaterial> material) { m_material = material; }

    /// Draw ground that edits moved by more than `threshold` (m) with another material, such as a darker one for
    /// freshly exposed regolith. Pass null to draw all ground with the ground material. Applies to bricks meshed
    /// after the call, so set it before the first update.
    void SetDisturbedMaterial(std::shared_ptr<ChVisualMaterial> material, double threshold = 0.005);

    /// Give the ground collision meshes of the given contact material, swept outward by `thickness` (m), in
    /// collision family `family` (0 to 15). Call before the first update.
    void EnableCollision(std::shared_ptr<ChContactMaterial> material, int family = 15, double thickness = 0.005);

    bool IsCollisionEnabled() const { return m_contact_material != nullptr; }
    int GetCollisionFamily() const { return m_family; }

    /// Keep the visual shapes in the visual model of GetBody() (default: true), where renderers that read the
    /// system's visual models, such as Chrono::Sensor, find them. Turn it off when a plugin that follows the
    /// shapes as they change, such as ChSiteVolumeVisualizationVSG, draws them instead: the VSG visual system binds
    /// a body's shapes once and would keep drawing the first ones.
    void SetUseVisualModel(bool val) { m_use_visual_model = val; }

    /// Mesh the bricks the volume changed since the last update and replace their shapes. Returns the number of
    /// bricks rebuilt.
    size_t Update();

    /// The fixed body holding the visual shapes.
    std::shared_ptr<ChBody> GetBody() const { return m_body; }
    std::shared_ptr<ChSiteVolume> GetVolume() const { return m_volume; }
    size_t GetNumShapes() const { return m_bricks.size(); }
    /// The visual shape of every brick holding ground; a rebuilt brick gets a new shape.
    std::vector<std::shared_ptr<ChVisualShapeTriangleMesh>> GetVisualShapes() const;
    /// Incremented each time an update replaces shapes.
    std::uint64_t GetShapesVersion() const { return m_shapes_version; }
    size_t GetNumTriangles() const { return m_num_triangles; }

  private:
    struct BrickShapes {
        std::shared_ptr<ChVisualShapeTriangleMesh> visual;
        std::shared_ptr<ChBody> collision;
        size_t triangles = 0;
    };

    void RemoveCollision(BrickShapes& brick);

    ChSystem* m_system;
    std::shared_ptr<ChSiteVolume> m_volume;
    std::shared_ptr<ChBody> m_body;
    std::shared_ptr<ChVisualMaterial> m_material;
    std::shared_ptr<ChVisualMaterial> m_disturbed_material;
    double m_disturbed_threshold = 0.005;
    std::shared_ptr<ChContactMaterial> m_contact_material;
    int m_family = 15;
    double m_thickness = 0.005;
    std::unordered_map<std::int64_t, BrickShapes> m_bricks;  ///< by brick
    std::uint64_t m_version = 0;                             ///< volume version meshed so far
    std::uint64_t m_shapes_version = 0;
    bool m_use_visual_model = true;
    size_t m_num_triangles = 0;
};

/// @} planet_module

}  // namespace planet
}  // namespace chrono

#endif
