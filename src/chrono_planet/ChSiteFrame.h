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
// East-north-up tangent frame at a landing site, mapping planet longitude,
// latitude and elevation to local Cartesian meters.
//
// =============================================================================

#ifndef CH_SITE_FRAME_H
#define CH_SITE_FRAME_H

#include <cmath>

#include "chrono/core/ChVector3.h"

#include "chrono_planet/ChApiPlanet.h"
#include "chrono_planet/core/Planet.h"

namespace chrono {
namespace planet {

/// @addtogroup planet_module
/// @{

/// Local east-north-up frame tangent to the planet sphere at a site origin.
/// Positions are equirectangular meters about the origin, so the frame is only meant for the
/// few kilometers a surface simulation covers; z is elevation above the origin's elevation.
class CH_PLANET_API ChSiteFrame {
  public:
    /// Planet radius in meters, shared with the terrain core.
    static constexpr double kRadius = qtplanet::kRadiusM;

    /// Construct a site frame at the given origin (degrees, meters).
    ChSiteFrame(double origin_lon_deg, double origin_lat_deg, double origin_elev_m)
        : m_lon0(origin_lon_deg), m_lat0(origin_lat_deg), m_elev0(origin_elev_m) {}

    /// Convert a planet position (degrees, meters) to site coordinates.
    ChVector3d ToLocal(double lon_deg, double lat_deg, double elev_m) const {
        const double y = Deg2Rad(lat_deg - m_lat0) * kRadius;
        const double x = Deg2Rad(WrapDelta(lon_deg - m_lon0)) * kRadius * std::cos(Deg2Rad(lat_deg));
        return ChVector3d(x, y, elev_m - m_elev0);
    }

    /// Convert site x/y meters to longitude and latitude in degrees.
    void ToLonLat(double x, double y, double& lon_deg, double& lat_deg) const {
        lat_deg = m_lat0 + Rad2Deg(y / kRadius);
        lon_deg = m_lon0 + Rad2Deg(x / (kRadius * std::cos(Deg2Rad(lat_deg))));
    }

    double GetOriginLongitude() const { return m_lon0; }  ///< origin longitude (degrees)
    double GetOriginLatitude() const { return m_lat0; }   ///< origin latitude (degrees)
    double GetOriginElevation() const { return m_elev0; } ///< origin elevation (meters)

  private:
    static double Deg2Rad(double d) { return d * (qtplanet::kPi / 180.0); }
    static double Rad2Deg(double r) { return r * (180.0 / qtplanet::kPi); }
    static double WrapDelta(double d) {
        while (d > 180.0)
            d -= 360.0;
        while (d < -180.0)
            d += 360.0;
        return d;
    }

    double m_lon0;
    double m_lat0;
    double m_elev0;
};

/// @} planet_module

}  // namespace planet
}  // namespace chrono

#endif
