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
// Solids given by signed distance functions, for editing a ChSiteVolume: the
// cutting geometry of an excavation tool, or a shape to add or remove.
//
// =============================================================================

#ifndef CH_SDF_SHAPE_H
#define CH_SDF_SHAPE_H

#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "chrono/core/ChFrame.h"
#include "chrono/core/ChVector3.h"
#include "chrono/geometry/ChAABB.h"

#include "chrono_planet/ChApiPlanet.h"

namespace chrono {
namespace planet {

/// @addtogroup planet_module
/// @{

/// A solid given by its signed distance function in its own frame: negative inside, positive outside. The
/// distance must be exact near the surface; farther away it may underestimate, as a union's does, since a
/// ChSiteVolume keeps distances only in a narrow band around its surface.
class CH_PLANET_API ChSdfShape {
  public:
    virtual ~ChSdfShape();

    /// Signed distance (m) from a point in the shape's frame.
    virtual double Distance(const ChVector3d& p) const = 0;

    /// Box holding the solid, in the shape's frame.
    virtual ChAABB GetBoundingBox() const = 0;
};

/// Sphere of the given radius, centered at the origin.
class CH_PLANET_API ChSdfSphere : public ChSdfShape {
  public:
    explicit ChSdfSphere(double radius);
    virtual double Distance(const ChVector3d& p) const override;
    virtual ChAABB GetBoundingBox() const override;

  private:
    double m_radius;
};

/// Box of the given side lengths centered at the origin, its edges rounded to `rounding` (m).
class CH_PLANET_API ChSdfBox : public ChSdfShape {
  public:
    explicit ChSdfBox(const ChVector3d& lengths, double rounding = 0);
    virtual double Distance(const ChVector3d& p) const override;
    virtual ChAABB GetBoundingBox() const override;

  private:
    ChVector3d m_half;  // half lengths less the rounding
    double m_rounding;
};

/// Capsule along z: a segment from -length/2 to length/2, swept by a sphere of the given radius.
class CH_PLANET_API ChSdfCapsule : public ChSdfShape {
  public:
    ChSdfCapsule(double radius, double length);
    virtual double Distance(const ChVector3d& p) const override;
    virtual ChAABB GetBoundingBox() const override;

  private:
    double m_radius;
    double m_half;
};

/// Cylinder along z of the given radius and length, centered at the origin.
class CH_PLANET_API ChSdfCylinder : public ChSdfShape {
  public:
    ChSdfCylinder(double radius, double length);
    virtual double Distance(const ChVector3d& p) const override;
    virtual ChAABB GetBoundingBox() const override;

  private:
    double m_radius;
    double m_half;
};

/// Union of shapes, each placed by a frame in the union's frame.
class CH_PLANET_API ChSdfUnion : public ChSdfShape {
  public:
    ChSdfUnion() = default;
    void AddShape(std::shared_ptr<ChSdfShape> shape, const ChFrame<>& frame = ChFrame<>());
    virtual double Distance(const ChVector3d& p) const override;
    virtual ChAABB GetBoundingBox() const override;

  private:
    std::vector<std::pair<std::shared_ptr<ChSdfShape>, ChFrame<>>> m_shapes;
};

/// A shape from a user function and its bounding box.
class CH_PLANET_API ChSdfFunction : public ChSdfShape {
  public:
    using Function = std::function<double(const ChVector3d& p)>;
    ChSdfFunction(Function distance, const ChAABB& box);
    virtual double Distance(const ChVector3d& p) const override;
    virtual ChAABB GetBoundingBox() const override;

  private:
    Function m_distance;
    ChAABB m_box;
};

/// @} planet_module

}  // namespace planet
}  // namespace chrono

#endif
