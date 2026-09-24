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

#include "chrono_planet/volume/ChSiteVolumeShapes.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "chrono/collision/ChCollisionShapeTriangleMesh.h"
#include "chrono/collision/ChCollisionSystem.h"

namespace chrono {
namespace planet {

namespace {
std::int64_t Key(const ChVector3i& b) {
    return (std::int64_t(b.x()) << 42) | (std::int64_t(b.y()) << 21) | std::int64_t(b.z());
}
}  // namespace

ChSiteVolumeShapes::ChSiteVolumeShapes(ChSystem* sys, std::shared_ptr<ChSiteVolume> volume) : m_system(sys), m_volume(std::move(volume)) {
    if (!m_system || !m_volume)
        throw std::invalid_argument("ChSiteVolumeShapes: null system or volume");
    m_body = chrono_types::make_shared<ChBody>();
    m_body->SetFixed(true);
    m_body->EnableCollision(false);
    m_system->AddBody(m_body);
    m_material = chrono_types::make_shared<ChVisualMaterial>();
    m_material->SetDiffuseColor(ChColor(0.5f, 0.5f, 0.5f));
}

ChSiteVolumeShapes::~ChSiteVolumeShapes() {
    for (auto& entry : m_bricks)
        RemoveCollision(entry.second);
}

void ChSiteVolumeShapes::SetDisturbedMaterial(std::shared_ptr<ChVisualMaterial> material, double threshold) {
    m_disturbed_material = material;
    m_disturbed_threshold = threshold;
}

void ChSiteVolumeShapes::EnableCollision(std::shared_ptr<ChContactMaterial> material, int family, double thickness) {
    if (family < 0 || family > 15)
        throw std::invalid_argument("ChSiteVolumeShapes: collision family must be in [0, 15]");
    m_contact_material = material;
    m_family = family;
    m_thickness = thickness;
}

void ChSiteVolumeShapes::RemoveCollision(BrickShapes& brick) {
    if (brick.collision) {
        m_system->RemoveBody(brick.collision);
        brick.collision.reset();
    }
}

size_t ChSiteVolumeShapes::Update() {
    if (m_volume->GetVersion() == m_version)
        return 0;
    std::vector<ChVector3i> changed;
    m_volume->GetChangedBricks(m_version, changed);
    m_version = m_volume->GetVersion();

    for (const auto& b : changed) {
        const std::int64_t key = Key(b);
        auto mesh = chrono_types::make_shared<ChTriangleMeshConnected>();
        const bool has_ground = m_volume->MeshBrick(b, *mesh, m_disturbed_threshold);
        auto it = m_bricks.find(key);
        if (it != m_bricks.end()) {
            m_num_triangles -= it->second.triangles;
            RemoveCollision(it->second);
            if (!has_ground) {
                m_bricks.erase(it);
                continue;
            }
        } else if (!has_ground) {
            continue;
        }
        BrickShapes& shapes = m_bricks[key];

        // A new shape for the new mesh: renderers key what they upload by the mesh, so an edited brick is rebuilt
        // and the others are left alone
        auto visual = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
        const auto& materials = mesh->GetIndicesMaterials();
        const bool disturbed = m_disturbed_material && std::any_of(materials.begin(), materials.end(), [](int m) { return m == 1; });
        visual->SetMesh(mesh);
        visual->SetMutable(false);
        visual->AddMaterial(m_material);
        if (disturbed)
            visual->AddMaterial(m_disturbed_material);
        else
            mesh->GetIndicesMaterials().clear();
        shapes.visual = visual;
        shapes.triangles = mesh->GetNumTriangles();
        m_num_triangles += shapes.triangles;

        if (m_contact_material) {
            auto body = chrono_types::make_shared<ChBody>();
            body->SetFixed(true);
            body->AddCollisionShape(chrono_types::make_shared<ChCollisionShapeTriangleMesh>(m_contact_material, mesh, true, false, m_thickness));
            body->EnableCollision(true);
            body->GetCollisionModel()->SetFamily(m_family);
            body->GetCollisionModel()->DisallowCollisionsWith(m_family);
            m_system->AddBody(body);
            // Collision models are bound at the system's first step; a body added later must be bound here
            if (auto* coll = m_system->GetCollisionSystem().get())
                coll->BindItem(body);
            shapes.collision = body;
        }
    }

    ++m_shapes_version;
    if (m_use_visual_model) {
        auto model = m_body->GetVisualModel();
        if (!model) {
            model = chrono_types::make_shared<ChVisualModel>();
            m_body->AddVisualModel(model);
        }
        model->Clear();
        for (const auto& entry : m_bricks)
            model->AddShape(entry.second.visual);
    }
    return changed.size();
}

std::vector<std::shared_ptr<ChVisualShapeTriangleMesh>> ChSiteVolumeShapes::GetVisualShapes() const {
    std::vector<std::shared_ptr<ChVisualShapeTriangleMesh>> shapes;
    shapes.reserve(m_bricks.size());
    for (const auto& entry : m_bricks)
        shapes.push_back(entry.second.visual);
    return shapes;
}

}  // namespace planet
}  // namespace chrono
