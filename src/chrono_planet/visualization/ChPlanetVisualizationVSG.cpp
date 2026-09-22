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

#include <algorithm>

#include "chrono/geometry/ChTriangleMeshConnected.h"

#include "chrono_vsg/utils/ChShapeBuilderVSG.h"

#include "chrono_planet/core/SphereMath.h"
#include "chrono_planet/lod/TileMeshBuilder.h"
#include "chrono_planet/visualization/ChPlanetVisualizationVSG.h"

namespace chrono {
namespace planet {

ChPlanetVisualizationVSG::ChPlanetVisualizationVSG(std::shared_ptr<QuadtreeWorld> world, const ChSiteFrame& site)
    : m_world(std::move(world)), m_site(site), m_wireframe(false), m_visible(true), m_mesh_version(~0ull) {
    m_material = chrono_types::make_shared<ChVisualMaterial>();
    m_material->SetDiffuseColor(ChColor(0.42f, 0.41f, 0.39f));
    m_material->SetRoughness(0.95f);
    m_material->SetMetallic(0.0f);
}

ChPlanetVisualizationVSG::~ChPlanetVisualizationVSG() {}

void ChPlanetVisualizationVSG::OnAttach() {
    m_vsys->SetCameraVertical(CameraVerticalDir::Z);
}

void ChPlanetVisualizationVSG::OnBindAssets() {
    m_terrain_scene = vsg::Switch::create();
    m_vsys->GetVSGScene()->addChild(m_terrain_scene);
}

void ChPlanetVisualizationVSG::OnRender() {
    // The camera in site coordinates, mapped onto the sphere for the LOD.
    const ChVector3d cam = m_vsys->GetCameraPosition();
    double lon, lat;
    m_site.ToLonLat(cam.x(), cam.y(), lon, lat);
    const double radius = qtplanet::kRadiusM + m_site.GetOriginElevation() + cam.z();
    const qtplanet::Vec3 p = qtplanet::pointOnSphere(lon, lat, radius);
    m_world->update(p.x, p.y, p.z, m_vsys->GetSimulationTime());

    SyncTiles();
}

void ChPlanetVisualizationVSG::SetVisible(bool val) {
    m_visible = val;
    if (m_terrain_scene)
        m_terrain_scene->setAllChildren(val);
}

void ChPlanetVisualizationVSG::SyncTiles() {
    const auto version = m_world->meshSetVersion();
    if (version == m_mesh_version)
        return;
    m_mesh_version = version;

    const auto live = m_world->getAllMeshes();

    // Drop tiles whose node is gone or whose mesh was regenerated.
    for (auto it = m_tiles.begin(); it != m_tiles.end();) {
        const auto found = live.find(it->first);
        if (found == live.end() || found->second->id != it->second.mesh_id) {
            auto& children = m_terrain_scene->children;
            const auto node = it->second.node;
            children.erase(std::remove_if(children.begin(), children.end(),
                                          [&](const vsg::Switch::Child& c) { return c.node == node; }),
                           children.end());
            it = m_tiles.erase(it);
        } else {
            ++it;
        }
    }

    // Add tiles the world holds that the scene does not.
    for (const auto& entry : live) {
        if (m_tiles.count(entry.first))
            continue;
        auto node = BuildTile(*entry.second);
        m_terrain_scene->addChild(m_visible, node);
        m_tiles[entry.first] = Tile{entry.second->id, node};
    }
}

vsg::ref_ptr<vsg::Node> ChPlanetVisualizationVSG::BuildTile(const Mesh& mesh) const {
    const size_t n = mesh.vertexData.size() / Mesh::kFloatsPerVertex;

    auto trimesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    auto& vertices = trimesh->GetCoordsVertices();
    auto& normals = trimesh->GetCoordsNormals();
    auto& uvs = trimesh->GetCoordsUV();
    vertices.reserve(n);
    normals.reserve(n);
    uvs.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        const float* v = &mesh.vertexData[i * Mesh::kFloatsPerVertex];
        // Center-relative planet position, re-projected through the site frame so it matches the physics.
        const double px = mesh.centerX + v[0];
        const double py = mesh.centerY + v[1];
        const double pz = mesh.centerZ + v[2];
        const qtplanet::LonLat ll = qtplanet::lonLatOf(px, py, pz);
        vertices.push_back(m_site.ToLocal(ll.lon, ll.lat, qtplanet::elevationOf(px, py, pz)));
        // Slopes are the normal's east and north components over its up component.
        ChVector3d normal(v[3], v[4], 1.0);
        normal.Normalize();
        normals.push_back(normal);
        uvs.push_back(ChVector2d(v[5], v[6]));
    }

    const auto& idx = MeshTopology::indices(mesh.level);
    auto& faces = trimesh->GetIndicesVertices();
    auto& face_normals = trimesh->GetIndicesNormals();
    auto& face_uvs = trimesh->GetIndicesUV();
    const size_t ntri = idx.size() / 3;
    faces.reserve(ntri);
    face_normals.reserve(ntri);
    face_uvs.reserve(ntri);
    for (size_t t = 0; t < ntri; ++t) {
        const ChVector3i tri(static_cast<int>(idx[3 * t]), static_cast<int>(idx[3 * t + 1]),
                             static_cast<int>(idx[3 * t + 2]));
        faces.push_back(tri);
        face_normals.push_back(tri);
        face_uvs.push_back(tri);
    }

    auto transform = vsg::MatrixTransform::create();
    return m_vsys->GetVSGShapeBuilder()->CreateTrimeshPbrMatShape(trimesh, transform, {m_material}, true, m_wireframe);
}

}  // namespace planet
}  // namespace chrono
