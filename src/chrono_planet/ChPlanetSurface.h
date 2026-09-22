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
// Deterministic planet surface: the DEM stack plus procedural base relief and
// craters, evaluated at one fixed level of detail.
//
// =============================================================================

#ifndef CH_PLANET_SURFACE_H
#define CH_PLANET_SURFACE_H

#include <memory>
#include <string>
#include <vector>

#include "chrono_planet/ChApiPlanet.h"

class MultiGeoTIFFManager;

namespace chrono {
namespace planet {

/// @addtogroup planet_module
/// @{

/// Height source for physics: DEM stack, base relief and craters at a fixed zoom.
/// The quadtree's own elevation depends on which tiles are resident, so a simulation samples
/// this instead; it evaluates the same terms the renderer draws at its finest level, with a
/// crater spacing matched to that level's mesh, so wheels and pixels ride one surface.
class CH_PLANET_API ChPlanetSurface {
  public:
    /// A GeoTIFF and the quadtree zoom range it answers for.
    struct DemSource {
        std::string path;  ///< absolute path, or relative to the working directory
        int min_zoom;      ///< coarsest zoom served by this DEM
        int max_zoom;      ///< finest zoom served by this DEM
    };

    /// Construct the surface over the given DEMs at a fixed zoom.
    /// A crater spacing of zero derives the spacing from the zoom, matching the renderer.
    /// A DEM that fails to open is logged and skipped; its area falls back to procedural noise.
    ChPlanetSurface(const std::vector<DemSource>& sources, int zoom, double crater_spacing_deg = 0.0);
    ~ChPlanetSurface();

    ChPlanetSurface(const ChPlanetSurface&) = delete;
    ChPlanetSurface& operator=(const ChPlanetSurface&) = delete;

    /// Surface elevation in meters at a longitude and latitude in degrees.
    double GetElevation(double lon_deg, double lat_deg) const;

    /// Fill out with an n x n row-major grid of elevations starting at (lon0, lat0) with the given steps.
    void GetElevationGrid(double lon0, double lat0, double step_lon, double step_lat, int n, std::vector<double>& out) const;

    int GetZoom() const { return m_zoom; }              ///< fixed sampling zoom
    double GetCraterSpacing() const;                    ///< crater spacing in degrees
    const MultiGeoTIFFManager& GetDemStack() const { return *m_dems; }  ///< underlying DEM stack

  private:
    double DemFallback(double lon_deg, double lat_deg) const;

    std::shared_ptr<MultiGeoTIFFManager> m_dems;
    int m_zoom;
    double m_crater_spacing_deg;
};

/// @} planet_module

}  // namespace planet
}  // namespace chrono

#endif
