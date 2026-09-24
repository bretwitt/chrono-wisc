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

#include "chrono_planet/visualization/ChSoilParticlesVisualizationVSG.h"

#include <stdexcept>

namespace chrono {
namespace planet {

ChSoilParticlesVisualizationVSG::ChSoilParticlesVisualizationVSG(std::shared_ptr<ChSoilParticles> particles) : m_particles(std::move(particles)) {
    if (!m_particles)
        throw std::invalid_argument("ChSoilParticlesVisualizationVSG: null particles");
}

ChSoilParticlesVisualizationVSG::~ChSoilParticlesVisualizationVSG() {}

void ChSoilParticlesVisualizationVSG::OnBindAssets() {
    m_scene = vsg::Switch::create();
    m_vsys->GetVSGScene()->addChild(m_scene);
    m_draws = std::make_unique<ChMeshDrawsVSG>(m_vsys, m_scene);
}

void ChSoilParticlesVisualizationVSG::SetVisible(bool val) {
    m_visible = val;
    if (m_draws)
        m_draws->SetVisible(val);
    if (m_scene)
        m_scene->setAllChildren(val);
}

void ChSoilParticlesVisualizationVSG::OnRender() {
    if (!m_draws || m_particles->GetMeshVersion() == m_version)
        return;
    m_version = m_particles->GetMeshVersion();
    m_draws->Remove(m_current);
    m_current.clear();
    const auto mesh = m_particles->GetMesh();
    if (mesh && mesh->GetNumTriangles() > 0)
        m_current = m_draws->Add(*mesh, {m_particles->GetMaterial()});
}

}  // namespace planet
}  // namespace chrono
