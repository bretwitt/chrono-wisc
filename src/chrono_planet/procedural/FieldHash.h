#ifndef QTPLANET_FIELDHASH_H
#define QTPLANET_FIELDHASH_H

#include "chrono_planet/ChApiPlanet.h"

// Grid-cell hash and octave fade shared by the procedural fields. Inline because they run per sample per octave.

#include <algorithm>

#include "chrono_planet/core/MathUtil.h"
#include "chrono_planet/core/Planet.h"

namespace qtfield {

// Ground kilometers per degree, the unit the fields' cell sizes are in.
constexpr double kKmPerDeg = qtplanet::kMetresPerDegree / 1000.0;

// Deterministic 32-bit hash of a grid cell and a salt.
inline unsigned hash32(int x, int y, int salt) {
    unsigned h = static_cast<unsigned>(x) * 374761393u ^ static_cast<unsigned>(y) * 668265263u ^
                 static_cast<unsigned>(salt) * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}

// hash32 mapped onto [0, 1).
inline float hash01(int x, int y, int salt) { return hash32(x, y, salt) / 4294967296.f; }

// floor() for doubles within int range.
inline int fastFloor(double x) {
    const int i = static_cast<int>(x);
    return i - (x < i ? 1 : 0);
}

// Returns an octave weight in [0, 1] as its feature becomes resolvable.
// fullAtSamples must exceed zeroBelowSamples; both thresholds are sample counts.
inline double resolutionFade(double samplesAcrossFeature, double zeroBelowSamples, double fullAtSamples) {
    if (samplesAcrossFeature < zeroBelowSamples) {
        return 0.0;
    }
    return qtplanet::smoothstep01(std::min(1.0, (samplesAcrossFeature - zeroBelowSamples) / (fullAtSamples - zeroBelowSamples)));
}

// C1 bump on |t| < 1, peak 1 at t = 0.
inline double bump(double t) {
    const double u = 1.0 - t * t;
    return (u > 0.0) ? u * u : 0.0;
}

}   // namespace qtfield

#endif   // QTPLANET_FIELDHASH_H
