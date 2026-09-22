#ifndef QTPLANET_FIELDGRID_H
#define QTPLANET_FIELDGRID_H

#include "chrono_planet/ChApiPlanet.h"

namespace qtfield {

// Square geographic sampling grid. Angles are degrees; steps must be positive
// and the grid must not cross the dateline. samplesPerSide must be positive.
// Sample (column, row) is (originLonDeg + column * stepLonDeg,
//                          originLatDeg + row * stepLatDeg).
struct CH_PLANET_API GridSpec {
    double originLonDeg;
    double originLatDeg;
    double stepLonDeg;
    double stepLatDeg;
    int samplesPerSide;
    // Positive angular resolution used to filter field octaves. This can differ
    // from the grid steps, e.g. when matching a coarser terrain LOD.
    double sampleSpacingDeg;
};

// Half-open row range: 0 <= begin <= end <= grid.samplesPerSide.
struct CH_PLANET_API RowRange {
    int begin;
    int end;
};

}   // namespace qtfield

#endif   // QTPLANET_FIELDGRID_H
