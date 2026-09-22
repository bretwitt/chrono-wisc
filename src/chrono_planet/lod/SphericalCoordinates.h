#ifndef QTPLANET_COORDINATE_TRAITS_SPHERICAL_H
#define QTPLANET_COORDINATE_TRAITS_SPHERICAL_H

#include "chrono_planet/ChApiPlanet.h"

#include <algorithm>
#include <array>
#include <memory>
#include <utility>
#include <vector>

#include "chrono_planet/lod/CoordinateSystems.h"
#include "chrono_planet/core/MathUtil.h"
#include "chrono_planet/core/Planet.h"
#include "chrono_planet/core/SphereMath.h"
#include "chrono_planet/lod/TileMetadata.h"

class MultiGeoTIFFManager;

inline constexpr bool kProceduralCraters = true;   // adds CraterField relief on top of the DEM
inline constexpr double kRootTileDeg = 16.0;       // QuadtreeWorld tile size, for spacing estimates

// Traits for the lon/lat sphere, covering tile bounds, relief sampling and cartesian grids.
template <>
class CoordinateTraits<Spherical> {
public:
    using Boundary = Spherical::Boundary;
    using Position = Spherical::Position;

    // NE, NW, SE, SW.
    static std::array<Boundary, 4> getChildBounds(const Boundary& b);

    // Chord distance from a camera position (meters) to the nearest point of the tile's elevation shell.
    [[nodiscard]] static double distanceToBounds(const Boundary& b, const qtplanet::Vec3& cameraM,
                                                 double minElev, double maxElev);

    // DEM plus procedural relief at one point, at zoomLevel's crater spacing. Null geoLoader means fallback terrain.
    [[nodiscard]] static double computeBaseElevation(const Position& pos, const MultiGeoTIFFManager* geoLoader,
                                                     int zoomLevel = 0);
    // CraterField height at a point for craters spaced spacingDeg, 0 when craters are off.
    [[nodiscard]] static double proceduralCraterRelief(double lonDeg, double latDeg, double spacingDeg);
    // The low-frequency Perlin swell under everything, meters.
    [[nodiscard]] static double proceduralBaseRelief(double lonDeg, double latDeg);

    // Cartesian positions of a (divisions+1)^2 grid over b, row-major south to north.
    // craterSpacingDeg <= 0 uses the grid's own spacing.
    [[nodiscard]] static std::vector<double> cartesianGrid(const Boundary& b, int divisions, const MultiGeoTIFFManager* geoLoader,
                                                           int zoomLevel, double craterSpacingDeg = 0.0);

    // Root tile grid indices (TileKey x, y) of the tile containing pos.
    [[nodiscard]] static TileKey computeTileIndices(const Position& pos, double tileSizeDegrees);
    // (lon, lat) of a root tile's center, degrees.
    [[nodiscard]] static qtplanet::LonLat tileCenterPosition(const TileKey& key, double tileSizeDegrees);

    static double wrapLongitude(double lon) { return qtplanet::wrapLongitude(lon); }
    // Height of a cartesian point above the sphere, meters.
    static double elevationOf(double x, double y, double z) { return qtplanet::elevationOf(x, y, z); }

    // Lowers a border vertex by depth meters toward the planet center.
    static qtplanet::Vec3 skirtVertex(const qtplanet::Vec3& p, double depth) {
        const double r = qtplanet::length(p);
        if (r < 1e-6) {
            return p;
        }
        const double s = (r - depth) / r;
        return {p.x * s, p.y * s, p.z * s};
    }

    // Clips a boundary to +-90 latitude and wraps its center longitude.
    static Boundary clampAtPoles(Boundary b) {
        b.centerLatDeg = std::clamp(b.centerLatDeg, -90.0, 90.0);
        if (b.centerLatDeg + b.halfHeightDeg > 90.0) {
            b.halfHeightDeg = 90.0 - b.centerLatDeg;
        }
        if (b.centerLatDeg - b.halfHeightDeg < -90.0) {
            b.halfHeightDeg = b.centerLatDeg + 90.0;
        }
        b.centerLonDeg = wrapLongitude(b.centerLonDeg);
        return b;
    }
};

#endif   // QTPLANET_COORDINATE_TRAITS_SPHERICAL_H
