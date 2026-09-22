#ifndef QTPLANET_SPHEREMATH_H
#define QTPLANET_SPHEREMATH_H

#include "chrono_planet/ChApiPlanet.h"

#include <cmath>

#include "chrono_planet/core/MathUtil.h"
#include "chrono_planet/core/Planet.h"

// Geometry of the body. Lon/lat to and from cartesian, the local horizon frame, and ground distances.

namespace qtplanet {

// Longitude folded into [-180, 180).
double wrapLongitude(double lonDeg);

// Floor on cos(lat) for anything walking a lon/lat grid, capping the polar stretch at 20x.
inline constexpr double kMinCosLat = 0.05;
// cos(lat) floored at kMinCosLat.
double cosLatClamped(double latDeg);

// Ground meters per degree. Not pole-guarded, so the longitude figure is exactly 0 at the poles.
constexpr double metresPerDegLat(double radiusM) { return radiusM * kPi / 180.0; }
double metresPerDegLon(double radiusM, double latDeg);

// Distance from the body's center.
inline double radiusOf(double x, double y, double z) { return std::sqrt(x * x + y * y + z * z); }
// Height above the reference sphere.
inline double elevationOf(double x, double y, double z) { return radiusOf(x, y, z) - kRadiusM; }

// Geographic position in degrees.
struct CH_PLANET_API LonLat {
    double lon, lat;
};

// Sub-point of a cartesian position, correct at any altitude.
LonLat lonLatOf(double x, double y, double z);
// Inverse of lonLatOf.
Vec3 pointOnSphere(double lonDeg, double latDeg, double radius);

// Orthonormal right-handed east, north, up at a point, in planet-centered cartesian coordinates.
struct CH_PLANET_API EnuFrame {
    Vec3 east, north, up;
};

// Frame at a lon/lat. Ill-defined at the poles.
EnuFrame enuAt(double lonDeg, double latDeg);
// Frame from an outward direction, with +x as east on the polar axis.
EnuFrame enuAlong(const Vec3& upDir);
// Outward radial at a lon/lat, the Up column of enuAt.
Vec3 dirFromLonLat(double lonDeg, double latDeg);

}   // namespace qtplanet

#endif   // QTPLANET_SPHEREMATH_H
