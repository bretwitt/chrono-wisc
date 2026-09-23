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

#include "chrono/core/ChTypes.h"

#include "chrono_planet/planets/mars/ChMars.h"

namespace chrono {
namespace planet {
namespace mars {

ChPlanetBody Body() {
    return ChPlanetBody("Mars", kRadius, kGravity);
}

ChCraterLayer::Params CraterParams() {
    ChCraterLayer::Params p;
    p.diameters_km = {80.0, 40.0, 20.0, 10.0, 5.0, 2.5, 1.25, 0.6, 0.3, 0.15, 0.08, 0.04, 0.02, 0.01};
    p.densities = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.3};
    // The transition diameter scales roughly with 1/g: 15 km on the Moon times 1.62 / 3.721.
    p.complex_diameter_km = 6.5;
    p.simple_depth_ratio = 0.196;
    p.simple_rim_ratio = 0.036;
    p.complex_depth_coef = 1.044;
    p.complex_depth_exp = 0.301;
    p.complex_rim_coef = 0.236;
    p.complex_rim_exp = 0.399;
    p.floor_flattening_km = 30.0;
    return p;
}

ChRockLayer::Params RockParams() {
    ChRockLayer::Params p;
    p.coverage = 0.05;
    return p;
}

std::shared_ptr<ChFilterChain> FilterChain() {
    const ChPlanetBody body = Body();
    auto chain = chrono_types::make_shared<ChFilterChain>();
    chain->AddFilter(chrono_types::make_shared<ChCraterLayer>(body, CraterParams()));
    chain->AddFilter(chrono_types::make_shared<ChRockLayer>(body, RockParams()));
    return chain;
}

std::shared_ptr<ChPlanetSurface> CreateSurface(const std::vector<ChGeoTiffSource>& dems, int zoom) {
    auto surface = chrono_types::make_shared<ChPlanetSurface>(Body(), zoom);
    for (const auto& dem : dems)
        surface->AddGeoTiff(dem);
    surface->SetFilterChain(FilterChain());
    return surface;
}

}  // namespace mars
}  // namespace planet
}  // namespace chrono
