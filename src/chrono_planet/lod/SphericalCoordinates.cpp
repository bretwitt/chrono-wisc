#include "chrono_planet/lod/SphericalCoordinates.h"

#include <algorithm>
#include <cmath>

#include "chrono_planet/ChPlanetSurface.h"
#include "chrono_planet/core/BenchProfiler.h"
#include "chrono_planet/core/Parallel.h"

namespace chrono {
namespace planet {

std::array<Spherical::Boundary, 4> CoordinateTraits<Spherical>::getChildBounds(const Boundary& b) {
    const double hw = b.halfWidthDeg * 0.5, hh = b.halfHeightDeg * 0.5;
    return {{clampAtPoles({wrapLongitude(b.centerLonDeg + hw), b.centerLatDeg + hh, hw, hh}),
             clampAtPoles({wrapLongitude(b.centerLonDeg - hw), b.centerLatDeg + hh, hw, hh}),
             clampAtPoles({wrapLongitude(b.centerLonDeg + hw), b.centerLatDeg - hh, hw, hh}),
             clampAtPoles({wrapLongitude(b.centerLonDeg - hw), b.centerLatDeg - hh, hw, hh})}};
}

double CoordinateTraits<Spherical>::distanceToBounds(const Boundary& b, const util::Vec3& cameraM,
                                                     double minElev, double maxElev, double sphereRadiusM) {
    // Clamp lon/lat to the rectangle and radius to the elevation range, then take the chord.
    const util::LonLat cam = util::lonLatOf(cameraM.x, cameraM.y, cameraM.z);
    const double camR = util::radiusOf(cameraM.x, cameraM.y, cameraM.z);
    const double elev = std::clamp(camR - sphereRadiusM, std::min(minElev, maxElev), std::max(minElev, maxElev));
    const double lat = std::clamp(cam.lat, b.centerLatDeg - b.halfHeightDeg, b.centerLatDeg + b.halfHeightDeg);
    const double lon = wrapLongitude(b.centerLonDeg + std::clamp(wrapLongitude(cam.lon - b.centerLonDeg), -b.halfWidthDeg, b.halfWidthDeg));
    const util::Vec3 surf = util::pointOnSphere(lon, lat, sphereRadiusM + elev);
    return util::length(cameraM - surf);
}

std::vector<double> CoordinateTraits<Spherical>::elevationGrid(const Boundary& b, int divisions, const ChPlanetSurface& surface,
                                                               int zoomLevel, double reliefSpacingDeg,
                                                               std::vector<double>* staticHeights) {
    const int n = divisions + 1;
    const size_t count = static_cast<size_t>(n) * n;
    const double stepLon = 2.0 * b.halfWidthDeg / divisions, stepLat = 2.0 * b.halfHeightDeg / divisions;
    const double startLat = b.centerLatDeg - b.halfHeightDeg;

    // Heights on the surface. The DEM sees the unwrapped origin; the surface wraps it for the relief.
    const double spacingDeg = reliefSpacingDeg > 0.0 ? reliefSpacingDeg : stepLon;
    const ChGeoGrid grid{b.centerLonDeg - b.halfWidthDeg, startLat, stepLon, stepLat, n, spacingDeg};
    std::vector<double> heights;
    if (!staticHeights) {
        surface.GetElevationGrid(grid, zoomLevel, heights);
        return heights;
    }
    if (staticHeights->size() == count) {
        heights = *staticHeights;
    } else {
        surface.GetStaticElevationGrid(grid, zoomLevel, heights);
        *staticHeights = heights;
    }
    surface.ApplyDynamicFilters(grid, heights);

    return heights;
}

std::vector<double> CoordinateTraits<Spherical>::cartesianGrid(const Boundary& b, int divisions, const ChPlanetSurface& surface,
                                                               int zoomLevel, double reliefSpacingDeg,
                                                               std::vector<double>* staticHeights) {
    const int n = divisions + 1;
    const size_t count = static_cast<size_t>(n) * n;
    const double stepLon = 2.0 * b.halfWidthDeg / divisions, stepLat = 2.0 * b.halfHeightDeg / divisions;
    const double startLon = wrapLongitude(b.centerLonDeg - b.halfWidthDeg), startLat = b.centerLatDeg - b.halfHeightDeg;

    const auto elev = elevationGrid(b, divisions, surface, zoomLevel, reliefSpacingDeg, staticHeights);

    // Trig of the grid's longitudes and latitudes.
    std::vector<double> cosLon(n), sinLon(n), cosLat(n), sinLat(n);
    for (int i = 0; i < n; ++i) {
        const double lam = util::deg2rad(wrapLongitude(startLon + i * stepLon));
        cosLon[i] = std::cos(lam);
        sinLon[i] = std::sin(lam);
    }
    for (int j = 0; j < n; ++j) {
        const double phi = util::deg2rad(startLat + j * stepLat);
        cosLat[j] = std::cos(phi);
        sinLat[j] = std::sin(phi);
    }

    const double radiusM = surface.GetBody().GetRadius();
    const bool par = util::parallelGrid(count);
    std::vector<double> outPos(count * 3);
    {
        BENCH_SCOPE("grid_sphere");
#pragma omp parallel for num_threads(util::parallelThreads()) schedule(static) if (par)
        for (int j = 0; j < n; ++j) {
            const double cp = cosLat[j], sp = sinLat[j];
            for (int i = 0; i < n; ++i) {
                const size_t v = static_cast<size_t>(j) * n + i;
                const double r = radiusM + elev[v], rcp = r * cp;   // == pointOnSphere(lon_i, lat_j, r)
                outPos[v * 3] = rcp * cosLon[i];
                outPos[v * 3 + 1] = rcp * sinLon[i];
                outPos[v * 3 + 2] = r * sp;
            }
        }
    }
    return outPos;
}

TileKey CoordinateTraits<Spherical>::computeTileIndices(const Position& pos, double tileSizeDegrees) {
    return {static_cast<int>(std::floor(wrapLongitude(pos.lon) / tileSizeDegrees)),
            static_cast<int>(std::floor((pos.lat + 90.0) / tileSizeDegrees))};
}

util::LonLat CoordinateTraits<Spherical>::tileCenterPosition(const TileKey& key, double tileSizeDegrees) {
    return {wrapLongitude((key.x + 0.5) * tileSizeDegrees), (key.y + 0.5) * tileSizeDegrees - 90.0};
}

}  // namespace planet
}  // namespace chrono
