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

#include "chrono_planet/volume/ChSdfShape.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace chrono {
namespace planet {

ChSdfShape::~ChSdfShape() {}

// -----------------------------------------------------------------------------

ChSdfSphere::ChSdfSphere(double radius) : m_radius(radius) {}

double ChSdfSphere::Distance(const ChVector3d& p) const {
    return p.Length() - m_radius;
}

ChAABB ChSdfSphere::GetBoundingBox() const {
    return ChAABB(ChVector3d(-m_radius), ChVector3d(m_radius));
}

// -----------------------------------------------------------------------------

ChSdfBox::ChSdfBox(const ChVector3d& lengths, double rounding) : m_rounding(std::max(0.0, rounding)) {
    for (int a = 0; a < 3; ++a)
        m_half[a] = std::max(0.0, 0.5 * lengths[a] - m_rounding);
}

double ChSdfBox::Distance(const ChVector3d& p) const {
    const ChVector3d q(std::abs(p.x()) - m_half.x(), std::abs(p.y()) - m_half.y(), std::abs(p.z()) - m_half.z());
    const ChVector3d outside(std::max(q.x(), 0.0), std::max(q.y(), 0.0), std::max(q.z(), 0.0));
    return outside.Length() + std::min(std::max(q.x(), std::max(q.y(), q.z())), 0.0) - m_rounding;
}

ChAABB ChSdfBox::GetBoundingBox() const {
    const ChVector3d h = m_half + ChVector3d(m_rounding);
    return ChAABB(-h, h);
}

// -----------------------------------------------------------------------------

ChSdfCapsule::ChSdfCapsule(double radius, double length) : m_radius(radius), m_half(0.5 * length) {}

double ChSdfCapsule::Distance(const ChVector3d& p) const {
    const double z = std::clamp(p.z(), -m_half, m_half);
    return (p - ChVector3d(0, 0, z)).Length() - m_radius;
}

ChAABB ChSdfCapsule::GetBoundingBox() const {
    return ChAABB(ChVector3d(-m_radius, -m_radius, -m_half - m_radius), ChVector3d(m_radius, m_radius, m_half + m_radius));
}

// -----------------------------------------------------------------------------

ChSdfCylinder::ChSdfCylinder(double radius, double length) : m_radius(radius), m_half(0.5 * length) {}

double ChSdfCylinder::Distance(const ChVector3d& p) const {
    const double dr = std::sqrt(p.x() * p.x() + p.y() * p.y()) - m_radius;
    const double dz = std::abs(p.z()) - m_half;
    const double outside = std::sqrt(std::max(dr, 0.0) * std::max(dr, 0.0) + std::max(dz, 0.0) * std::max(dz, 0.0));
    return outside + std::min(std::max(dr, dz), 0.0);
}

ChAABB ChSdfCylinder::GetBoundingBox() const {
    return ChAABB(ChVector3d(-m_radius, -m_radius, -m_half), ChVector3d(m_radius, m_radius, m_half));
}

// -----------------------------------------------------------------------------

void ChSdfUnion::AddShape(std::shared_ptr<ChSdfShape> shape, const ChFrame<>& frame) {
    m_shapes.push_back({std::move(shape), frame});
}

double ChSdfUnion::Distance(const ChVector3d& p) const {
    double d = std::numeric_limits<double>::max();
    for (const auto& [shape, frame] : m_shapes)
        d = std::min(d, shape->Distance(frame.TransformPointParentToLocal(p)));
    return d;
}

ChAABB ChSdfUnion::GetBoundingBox() const {
    ChAABB box;
    for (const auto& [shape, frame] : m_shapes)
        box += shape->GetBoundingBox().Transform(frame);
    return box;
}

// -----------------------------------------------------------------------------

ChSdfFunction::ChSdfFunction(Function distance, const ChAABB& box) : m_distance(std::move(distance)), m_box(box) {}

double ChSdfFunction::Distance(const ChVector3d& p) const {
    return m_distance(p);
}

ChAABB ChSdfFunction::GetBoundingBox() const {
    return m_box;
}

}  // namespace planet
}  // namespace chrono
