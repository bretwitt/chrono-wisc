#ifndef QTPLANET_COORDINATESYSTEMS_H
#define QTPLANET_COORDINATESYSTEMS_H

#include "chrono_planet/ChApiPlanet.h"

// The lon/lat coordinate system, with position and tile-bounds types in degrees.
struct CH_PLANET_API Spherical {
    // A point on the sphere, in degrees.
    struct Position {
        double lon, lat;
    };

    // Geographic center and half extents, all in degrees.
    struct Boundary {
        double centerLonDeg, centerLatDeg, halfWidthDeg, halfHeightDeg;
    };
};

// Geometry operations for one coordinate system, specialized in that system's own header.
template <typename CoordSystem>
struct CoordinateTraits;

#endif   // QTPLANET_COORDINATESYSTEMS_H
