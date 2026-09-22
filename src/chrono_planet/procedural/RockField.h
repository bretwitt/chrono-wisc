#ifndef QTPLANET_ROCKFIELD_H
#define QTPLANET_ROCKFIELD_H

#include "chrono_planet/ChApiPlanet.h"

#include <array>
#include <cstdint>
#include <vector>

#include "chrono_planet/procedural/FieldGrid.h"

// One boulder placement, deterministic from its grid cell.
struct CH_PLANET_API RockInstance {
    std::uint64_t id = 0;                   // stable and unique per rock
    double lonDeg = 0.0, latDeg = 0.0;      // center
    float radiusM = 0.f;                    // half the longest axis, meters
    float yawRad = 0.f;                     // spin about the local up
    float tiltRad = 0.f, tiltAzRad = 0.f;   // lean off vertical, and its azimuth
    float buryFrac = 0.f;                   // fraction of the diameter below ground
    std::uint16_t meshId = 0;               // index into RockMeshes
    std::uint8_t sizeClass = 0;             // 0 = finest class
};

// Golombek-Rapp boulder field, a pure function of (lon, lat).
class CH_PLANET_API RockField {
public:
    static constexpr int kClasses = 6;   // one octave each, finest first
    // Top-edge diameter (meters) of class sizeClass, 0 outside the range.
    [[nodiscard]] static double classDiameterM(int sizeClass);
    // Largest radius any rock can have, half the coarsest class's top edge.
    [[nodiscard]] static double maxRadiusM();
    // Rocks with radius >= minRadiusM in a lon/lat rect, half-open on the max edges.
    [[nodiscard]] static std::vector<RockInstance> query(double minLonDeg, double minLatDeg, double maxLonDeg, double maxLatDeg,
                                                         double minRadiusM);
    // Unit quaternion (x, y, z, w) for the rock's yaw and lean, in its local ENU frame.
    [[nodiscard]] static std::array<double, 4> orientation(const RockInstance& rock);
    // Height of the mesh origin above ground once buried per buryFrac, with a floor on the depth.
    [[nodiscard]] static double centreRiseM(const RockInstance& rock, double bottomExtentM,
                                            double totalHeightM);
    // Plinth relief (meters) at one point from every rock whose bedding reaches it.
    [[nodiscard]] static double beddingAt(double lonDeg, double latDeg, double sampleSpacingDeg);
    // Accumulates beddingAt() into heightsM only within rows [begin, end).
    // heightsM must contain at least samplesPerSide squared elements, row-major;
    // it is never resized. Concurrent calls must own disjoint row ranges.
    static void addBeddingGridRows(const qtfield::GridSpec& grid, qtfield::RowRange rows, std::vector<double>& heightsM);
};

#endif   // QTPLANET_ROCKFIELD_H
