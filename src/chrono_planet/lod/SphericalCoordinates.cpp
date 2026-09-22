#include "chrono_planet/lod/SphericalCoordinates.h"

#include <algorithm>
#include <cmath>

#include "chrono_planet/core/BenchProfiler.h"
#include "chrono_planet/core/Parallel.h"
#include "chrono_planet/procedural/CraterField.h"
#include "chrono_planet/dem/MultiGeoTIFFManager.h"
#include "chrono_planet/procedural/Perlin.h"
#include "chrono_planet/lod/TileMeshBuilder.h"

using qtplanet::kRadiusM;

namespace {

// Stand-in height where no DEM covers a point.
double fallbackElevation(double lonDeg, double latDeg) {
    return Perlin::noise(static_cast<float>(latDeg * 0.1), static_cast<float>(lonDeg * 0.1));
}

// The ~545 km swell under everything, low-frequency Perlin on the sphere.
constexpr float kBaseFreq = 20.0f;
constexpr double kBaseAmplitudeM = 1000.0;
// Tuned gain on the swell. Kept as a separate factor so the sum rounds as it always has.
constexpr double kBaseReliefGain = 0.7;
double baseReliefM(float sphereNoise) { return kBaseReliefGain * (sphereNoise * kBaseAmplitudeM); }

}   // namespace

std::array<Spherical::Boundary, 4> CoordinateTraits<Spherical>::getChildBounds(const Boundary& b) {
    const double hw = b.halfWidthDeg * 0.5, hh = b.halfHeightDeg * 0.5;
    return {{clampAtPoles({wrapLongitude(b.centerLonDeg + hw), b.centerLatDeg + hh, hw, hh}),
             clampAtPoles({wrapLongitude(b.centerLonDeg - hw), b.centerLatDeg + hh, hw, hh}),
             clampAtPoles({wrapLongitude(b.centerLonDeg + hw), b.centerLatDeg - hh, hw, hh}),
             clampAtPoles({wrapLongitude(b.centerLonDeg - hw), b.centerLatDeg - hh, hw, hh})}};
}

double CoordinateTraits<Spherical>::distanceToBounds(const Boundary& b, const qtplanet::Vec3& cameraM,
                                                     double minElev, double maxElev) {
    // Clamp lon/lat to the rectangle and radius to the elevation range, then take the chord.
    const qtplanet::LonLat cam = qtplanet::lonLatOf(cameraM.x, cameraM.y, cameraM.z);
    const double camR = qtplanet::radiusOf(cameraM.x, cameraM.y, cameraM.z);
    const double elev = std::clamp(camR - kRadiusM, std::min(minElev, maxElev), std::max(minElev, maxElev));
    const double lat = std::clamp(cam.lat, b.centerLatDeg - b.halfHeightDeg, b.centerLatDeg + b.halfHeightDeg);
    const double lon = wrapLongitude(b.centerLonDeg + std::clamp(wrapLongitude(cam.lon - b.centerLonDeg), -b.halfWidthDeg, b.halfWidthDeg));
    const qtplanet::Vec3 surf = qtplanet::pointOnSphere(lon, lat, kRadiusM + elev);
    return qtplanet::length(cameraM - surf);
}

double CoordinateTraits<Spherical>::proceduralBaseRelief(double lonDeg, double latDeg) {
    return baseReliefM(Perlin::onSphere(qtplanet::dirFromLonLat(lonDeg, latDeg), kBaseFreq));
}

double CoordinateTraits<Spherical>::proceduralCraterRelief(double lonDeg, double latDeg, double spacingDeg) {
    return kProceduralCraters ? CraterField::heightAt(lonDeg, latDeg, spacingDeg) : 0.0;
}

std::vector<double> CoordinateTraits<Spherical>::cartesianGrid(const Boundary& b, int divisions,
                                                               const MultiGeoTIFFManager* geoLoader,
                                                               int zoomLevel, double craterSpacingDeg) {
    const int n = divisions + 1;
    const size_t count = static_cast<size_t>(n) * n;
    const double stepLon = 2.0 * b.halfWidthDeg / divisions, stepLat = 2.0 * b.halfHeightDeg / divisions;
    const double startLon = wrapLongitude(b.centerLonDeg - b.halfWidthDeg), startLat = b.centerLatDeg - b.halfHeightDeg;

    // Cache the lon/lat and their trig for the grid, so the DEM and Perlin can be sampled in parallel.
    std::vector<double> lons(count), lats(count);
    std::vector<double> cosLon(n), sinLon(n), cosLat(n), sinLat(n);
    for (int i = 0; i < n; ++i) {
        lons[i] = wrapLongitude(startLon + i * stepLon);
        const double lam = qtplanet::deg2rad(lons[i]);
        cosLon[i] = std::cos(lam);
        sinLon[i] = std::sin(lam);
    }
    for (int j = 0; j < n; ++j) {
        const double lat = startLat + j * stepLat;
        const double phi = qtplanet::deg2rad(lat);
        cosLat[j] = std::cos(phi);
        sinLat[j] = std::sin(phi);
        for (int i = 0; i < n; ++i) {
            lons[j * n + i] = lons[i];
            lats[j * n + i] = lat;
        }
    }

    // The grid passes fork the OpenMP team for bake-sized grids and stay serial for the mesh lattice.
    const bool par = qtplanet::parallelGrid(count);
    std::vector<double> elev(count);
    std::optional<MultiGeoTIFFManager::GridSamples> sampled;
    {
        BENCH_SCOPE("grid_dem");
        if (geoLoader) {
            sampled = geoLoader->sampleGrid(b.centerLonDeg - b.halfWidthDeg, startLat, stepLon, stepLat, n, n, zoomLevel);
        }
    }
    if (sampled) {
        elev = std::move(sampled->elevations);
        for (size_t v : sampled->missing) {
            elev[v] = fallbackElevation(lons[v], lats[v]);
        }

        {
            BENCH_SCOPE("grid_perlin");
            // Rows are independent. Each thread keeps its own row scratch.
#pragma omp parallel num_threads(qtplanet::parallelThreads()) if (par)
            {
                std::vector<qtplanet::Vec3> dirs(n);
                std::vector<float> relief(n), scratch;
#pragma omp for schedule(static)
                for (int j = 0; j < n; ++j) {
                    const double cp = cosLat[j], sp = sinLat[j];
                    for (int i = 0; i < n; ++i) {
                        dirs[i] = {cp * cosLon[i], cp * sinLon[i], sp};   // == dirFromLonLat(lon_i, lat_j)
                    }
                    Perlin::writeSphereNoiseRow(dirs.data(), n, kBaseFreq, relief.data(), scratch);
                    double* row = &elev[static_cast<size_t>(j) * n];
                    for (int i = 0; i < n; ++i) {
                        row[i] += baseReliefM(relief[i]);
                    }
                }
            }
        }

        if (kProceduralCraters) {
            BENCH_SCOPE("grid_craters");
            const double spacingDeg = craterSpacingDeg > 0.0 ? craterSpacingDeg : stepLon;
            // Guard against the dateline wrap.
            if (wrapLongitude(startLon + (n - 1) * stepLon) >= startLon) {
                CraterField::addGrid({startLon, startLat, stepLon, stepLat, n, spacingDeg}, elev);
            } else {
#pragma omp parallel for num_threads(qtplanet::parallelThreads()) schedule(static) if (par)
                for (size_t v = 0; v < count; ++v) {
                    elev[v] += CraterField::heightAt(lons[v], lats[v], spacingDeg);
                }
            }
        }
    } else {
        BENCH_SCOPE("grid_pointwise");
        // sampleGrid is const and guards its OGR transforms with a mutex.
#pragma omp parallel for num_threads(qtplanet::parallelThreads()) schedule(dynamic, 64) if (par)
        for (size_t v = 0; v < count; ++v) {
            elev[v] = computeBaseElevation({lons[v], lats[v]}, geoLoader, zoomLevel);
        }
    }

    std::vector<double> outPos(count * 3);
    {
        BENCH_SCOPE("grid_sphere");
#pragma omp parallel for num_threads(qtplanet::parallelThreads()) schedule(static) if (par)
        for (int j = 0; j < n; ++j) {
            const double cp = cosLat[j], sp = sinLat[j];
            for (int i = 0; i < n; ++i) {
                const size_t v = static_cast<size_t>(j) * n + i;
                const double r = kRadiusM + elev[v], rcp = r * cp;   // == pointOnSphere(lon_i, lat_j, r)
                outPos[v * 3] = rcp * cosLon[i];
                outPos[v * 3 + 1] = rcp * sinLon[i];
                outPos[v * 3 + 2] = r * sp;
            }
        }
    }
    return outPos;
}

double CoordinateTraits<Spherical>::computeBaseElevation(const Position& pos,
                                                         const MultiGeoTIFFManager* geoLoader,
                                                         int zoomLevel) {
    if (!geoLoader) {
        return fallbackElevation(pos.lon, pos.lat);
    }
    const double dem = geoLoader->sample(pos.lon, pos.lat, zoomLevel).value_or(fallbackElevation(pos.lon, pos.lat));
    const double spacingDeg = MeshTopology::vertexSpacingDeg(zoomLevel, kRootTileDeg);
    return dem + proceduralBaseRelief(pos.lon, pos.lat) + proceduralCraterRelief(pos.lon, pos.lat, spacingDeg);
}

TileKey CoordinateTraits<Spherical>::computeTileIndices(const Position& pos, double tileSizeDegrees) {
    return {static_cast<int>(std::floor(wrapLongitude(pos.lon) / tileSizeDegrees)),
            static_cast<int>(std::floor((pos.lat + 90.0) / tileSizeDegrees))};
}

qtplanet::LonLat CoordinateTraits<Spherical>::tileCenterPosition(const TileKey& key, double tileSizeDegrees) {
    return {wrapLongitude((key.x + 0.5) * tileSizeDegrees), (key.y + 0.5) * tileSizeDegrees - 90.0};
}
