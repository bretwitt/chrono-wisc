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

#ifndef CH_API_PLANET_H
#define CH_API_PLANET_H

#include "chrono/ChVersion.h"
#include "chrono/core/ChPlatform.h"

// When compiling this library, remember to define CH_API_COMPILE_PLANET
// (so that the symbols with 'CH_PLANET_API' in front of them will be
// marked as exported). Otherwise, just do not define it if you
// link the library to your code, and the symbols will be imported.

#if defined(CH_API_COMPILE_PLANET)
    #define CH_PLANET_API ChApiEXPORT
#else
    #define CH_PLANET_API ChApiIMPORT
#endif

/**
    @defgroup planet_module PLANET module
    @brief Planetary terrain from GeoTIFF elevation models

    This module provides a headless planetary terrain core: a stack of GeoTIFF digital
    elevation models sampled at any longitude, latitude and level of detail, procedural
    crater, rock and roughness fields layered on top of it, and a quadtree of tiles that
    streams around a moving viewpoint. Chrono::Vehicle builds rigid and SCM terrains on
    it (see PlanetTerrain and PlanetSCMTerrain).

    For additional information, see:
    - the [installation guide](@ref module_planet_installation)
    - the [tutorials](@ref tutorial_root)

    @{
        @defgroup planet_core Body constants, sphere math and profiling
        @defgroup planet_dem GeoTIFF loading and the multi-DEM elevation stack
        @defgroup planet_procedural Craters, rocks and surface roughness
        @defgroup planet_lod Quadtree, tiles and tile mesh construction
        @defgroup planet_shading Hapke reflectance parameters
    @}
*/

namespace chrono {

/// @addtogroup planet_module
/// @{

/// Namespace with classes for the PLANET module.
namespace planet {}

/// @}

}  // namespace chrono

#endif
