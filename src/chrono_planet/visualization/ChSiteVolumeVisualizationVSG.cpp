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

#include "chrono_planet/visualization/ChSiteVolumeVisualizationVSG.h"

#include <stdexcept>
#include <unordered_set>

namespace chrono {
namespace planet {

ChSiteVolumeVisualizationVSG::ChSiteVolumeVisualizationVSG(std::shared_ptr<ChSiteVolumeShapes> shapes) : m_shapes(std::move(shapes)) {
    if (!m_shapes)
        throw std::invalid_argument("ChSiteVolumeVisualizationVSG: null shapes");
}

ChSiteVolumeVisualizationVSG::~ChSiteVolumeVisualizationVSG() {}

void ChSiteVolumeVisualizationVSG::OnBindAssets() {
    // Bricks are added from the first frame on, as the terrain plugin adds its tiles: state built while the visual
    // system is being initialized would not be compiled
    m_scene = vsg::Switch::create();
    m_vsys->GetVSGScene()->addChild(m_scene);
    m_draws = std::make_unique<ChMeshDrawsVSG>(m_vsys, m_scene);
}

void ChSiteVolumeVisualizationVSG::OnRender() {
    Sync();
}

void ChSiteVolumeVisualizationVSG::SetVisible(bool val) {
    m_visible = val;
    if (m_draws)
        m_draws->SetVisible(val);
    if (m_scene)
        m_scene->setAllChildren(val);
}

void ChSiteVolumeVisualizationVSG::Sync() {
    if (!m_scene)
        return;
    if (m_wireframe_requested != m_wireframe) {
        m_wireframe = m_wireframe_requested;
        m_draws->Clear();
        m_draws->SetWireframe(m_wireframe);
        m_bricks.clear();
        m_version = ~0ull;
    }
    if (m_shapes->GetShapesVersion() == m_version)
        return;
    m_version = m_shapes->GetShapesVersion();

    // Shapes are replaced, not edited, when their brick is rebuilt, so they are matched by identity
    const auto shapes = m_shapes->GetVisualShapes();
    std::unordered_set<const ChVisualShapeTriangleMesh*> live;
    for (const auto& shape : shapes)
        live.insert(shape.get());
    for (auto it = m_bricks.begin(); it != m_bricks.end();) {
        if (live.count(it->first)) {
            ++it;
            continue;
        }
        m_draws->Remove(it->second.draws);
        it = m_bricks.erase(it);
    }
    for (const auto& shape : shapes) {
        if (m_bricks.count(shape.get()))
            continue;
        BrickNodes& nodes = m_bricks[shape.get()];
        nodes.shape = shape;
        nodes.draws = m_draws->Add(*shape->GetMesh(), shape->GetMaterials());
    }
}

}  // namespace planet
}  // namespace chrono
