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
// Mars preset: body constants and a starting set of relief parameters.
//
// =============================================================================

#ifndef CH_MARS_H
#define CH_MARS_H

#include <memory>
#include <vector>

#include "chrono_planet/ChApiPlanet.h"
#include "chrono_planet/ChPlanetBody.h"
#include "chrono_planet/ChPlanetSurface.h"
#include "chrono_planet/dem/ChGeoTiffStack.h"
#include "chrono_planet/filters/ChSurfaceFilter.h"
#include "chrono_planet/procedural/ChCraterLayer.h"
#include "chrono_planet/procedural/ChRockLayer.h"

namespace chrono {
namespace planet {

/// Mars, as a sphere of its mean radius.
/// The relief parameters are NOT calibrated: they reuse the lunar crater shapes with the
/// simple-to-complex transition scaled by 1/gravity, and a moderate rock cover. They are a
/// starting point to tune against site data, not a validated model.
namespace mars {

/// @addtogroup planet_module
/// @{

constexpr double kRadius = 3389500.0;  ///< mean radius (m)
constexpr double kGravity = 3.721;     ///< surface gravity (m/s^2)

/// Mars as a body: radius, gravity, and a spherical geographic system of that radius.
/// DEMs referenced to the areoid (MOLA MEGDR, HiRISE DTMs) are heights relative to a surface within
/// a few kilometers of this sphere; set ChGeoTiffSource::offset if a site needs a common datum.
CH_PLANET_API ChPlanetBody Body();

/// Craters, uncalibrated: lunar morphometry with a 6.5 km simple-to-complex transition.
CH_PLANET_API ChCraterLayer::Params CraterParams();

/// Boulders, uncalibrated: Golombek-Rapp with k = 0.05.
CH_PLANET_API ChRockLayer::Params RockParams();

/// The preset filter chain: the crater and rock layers. Each call returns a new chain.
CH_PLANET_API std::shared_ptr<ChFilterChain> FilterChain();

/// A Martian surface: the Mars body, the given DEMs, and the preset filter chain.
/// Throws std::runtime_error if a DEM cannot be loaded.
CH_PLANET_API std::shared_ptr<ChPlanetSurface> CreateSurface(const std::vector<ChGeoTiffSource>& dems = std::vector<ChGeoTiffSource>(),
                                                             int zoom = 15);

/// @} planet_module

}  // namespace mars
}  // namespace planet
}  // namespace chrono

#endif
