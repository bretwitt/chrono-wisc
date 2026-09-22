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

#include "chrono_planet/ChPlanetSurface.h"

#include "chrono_planet/dem/MultiGeoTIFFManager.h"
#include "chrono_planet/lod/SphericalCoordinates.h"
#include "chrono_planet/lod/TileMeshBuilder.h"
#include "chrono_planet/procedural/CraterField.h"
#include "chrono_planet/procedural/Perlin.h"

namespace chrono {
namespace planet {

ChPlanetSurface::ChPlanetSurface(const std::vector<DemSource>& sources, int zoom, double crater_spacing_deg)
    : m_dems(std::make_shared<MultiGeoTIFFManager>()), m_zoom(zoom), m_crater_spacing_deg(crater_spacing_deg) {
    for (const auto& s : sources)
        (void)m_dems->addSource(s.path, s.min_zoom, s.max_zoom);  // a missing DEM is logged by the loader
}

ChPlanetSurface::~ChPlanetSurface() = default;

// Same fallback the tile builder uses where no DEM answers, so physics and renderer agree off-DEM.
double ChPlanetSurface::DemFallback(double lon_deg, double lat_deg) const {
    return Perlin::noise(lat_deg * 0.1, lon_deg * 0.1);
}

double ChPlanetSurface::GetElevation(double lon_deg, double lat_deg) const {
    lon_deg = CoordinateTraits<Spherical>::wrapLongitude(lon_deg);
    const double dem = m_dems->sample(lon_deg, lat_deg, m_zoom).value_or(DemFallback(lon_deg, lat_deg));
    return dem + CoordinateTraits<Spherical>::proceduralBaseRelief(lon_deg, lat_deg) +
           CoordinateTraits<Spherical>::proceduralCraterRelief(lon_deg, lat_deg, GetCraterSpacing());
}

void ChPlanetSurface::GetElevationGrid(double lon0,
                                       double lat0,
                                       double step_lon,
                                       double step_lat,
                                       int n,
                                       std::vector<double>& out) const {
    const size_t count = static_cast<size_t>(n) * n;
    out.assign(count, 0.0);

    std::vector<size_t> missing;
    if (auto grid = m_dems->sampleGrid(lon0, lat0, step_lon, step_lat, n, n, m_zoom)) {
        out = std::move(grid->elevations);
        missing = std::move(grid->missing);
    } else {
        missing.resize(count);
        for (size_t v = 0; v < count; ++v)
            missing[v] = v;
        out.assign(count, 0.0);
    }
    for (size_t v : missing) {
        const double lon = lon0 + (v % n) * step_lon;
        const double lat = lat0 + (v / n) * step_lat;
        out[v] = DemFallback(lon, lat);
    }

    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            out[static_cast<size_t>(j) * n + i] +=
                CoordinateTraits<Spherical>::proceduralBaseRelief(lon0 + i * step_lon, lat0 + j * step_lat);

    CraterField::addGrid({lon0, lat0, step_lon, step_lat, n, GetCraterSpacing()}, out);
}

// The true vertex spacing of a level-zoom tile mesh, so physics resolves the same crater octaves as the drawn mesh.
double ChPlanetSurface::GetCraterSpacing() const {
    if (m_crater_spacing_deg > 0.0)
        return m_crater_spacing_deg;
    return MeshTopology::vertexSpacingDeg(m_zoom, kRootTileDeg);
}

}  // namespace planet
}  // namespace chrono
