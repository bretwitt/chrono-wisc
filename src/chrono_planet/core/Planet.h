#ifndef QTPLANET_PLANET_H
#define QTPLANET_PLANET_H

#include "chrono_planet/ChApiPlanet.h"

// The body being rendered. The terrain and rock fragment shaders hard-code kRadiusM and change with it.

namespace qtplanet {

// Radius and relief envelope of the body being rendered.
struct CH_PLANET_API Body {
    double radiusM;
    // Terrain never reaches deeper than this below the mean radius. An underestimate culls visible ground.
    double deepestRelief;
};

// The Moon, 1737.4 km mean radius with 20 km of relief margin below it.
inline constexpr Body kMoon{1737400.0, 20000.0};

// The body the whole project renders.
inline constexpr Body kBody = kMoon;

inline constexpr double kRadiusM = kBody.radiusM;
inline constexpr double kPi = 3.14159265358979323846;
inline constexpr double kTwoPi = 2.0 * kPi;
// Ground meters per degree along a meridian (or the equator).
inline constexpr double kMetresPerDegree = kRadiusM * kPi / 180.0;
inline constexpr double kCircumferenceM = 2.0 * kPi * kRadiusM;

// Far plane. Six radii keep the limb in frame from the surface out to the flyin demo's opening altitude.
inline constexpr double kFarPlaneM = 6.0 * kRadiusM;

}   // namespace qtplanet

#endif   // QTPLANET_PLANET_H
