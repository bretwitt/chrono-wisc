#ifndef QTPLANET_SURFACEROUGHNESS_H
#define QTPLANET_SURFACEROUGHNESS_H

#include "chrono_planet/ChApiPlanet.h"

#include <vector>

#include "chrono_planet/procedural/FieldGrid.h"

// Micro-surface roughness below the DEM's band, resolved to the sample spacing.
class CH_PLANET_API SurfaceRoughness {
public:
    // Roughness (meters) at one point, with the octaves sampleSpacingDeg resolves.
    [[nodiscard]] static double heightAt(double lonDeg, double latDeg, double sampleSpacingDeg);
    // Accumulates heightAt() into heightsM only within rows [begin, end).
    // heightsM must contain at least samplesPerSide squared elements, row-major;
    // it is never resized. Concurrent calls must own disjoint row ranges.
    static void addGridRows(const qtfield::GridSpec& grid, qtfield::RowRange rows, std::vector<double>& heightsM);
    [[nodiscard]] static double slopePerOctave();   // RMS slope of the coarsest octave
    [[nodiscard]] static double fineBoost();        // slope multiplier per octave toward the fine end
    // RMS slope (rise/run) of the octaves a sampling at sampleSpacingDeg resolves.
    [[nodiscard]] static double rmsSlopeAt(double sampleSpacingDeg);
};

#endif   // QTPLANET_SURFACEROUGHNESS_H
