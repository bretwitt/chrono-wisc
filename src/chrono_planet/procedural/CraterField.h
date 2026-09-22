#ifndef QTPLANET_CRATERFIELD_H
#define QTPLANET_CRATERFIELD_H

#include "chrono_planet/ChApiPlanet.h"

#include <vector>

#include "chrono_planet/procedural/FieldGrid.h"

// Procedural crater field with a lunar equilibrium size distribution, a pure function of (lon, lat).
class CH_PLANET_API CraterField {
public:
    // Combined crater, rock-bedding and surface-roughness relief in meters.
    // lonDeg/latDeg are degrees; positive sampleSpacingDeg filters unresolved octaves.
    [[nodiscard]] static double heightAt(double lonDeg, double latDeg, double sampleSpacingDeg);

    // Accumulates relief into heightsM, including rock bedding and surface roughness.
    // heightsM must contain at least samplesPerSide squared elements, row-major.
    // Existing elevations are preserved and added to; the buffer is never resized.
    static void addGrid(const qtfield::GridSpec& grid, std::vector<double>& heightsM);
};

#endif   // QTPLANET_CRATERFIELD_H
