#ifndef QTPLANET_ASSETPATHS_H
#define QTPLANET_ASSETPATHS_H

#include "chrono_planet/ChApiPlanet.h"

#include <string>

namespace qtplanet {

// The directory holding shaders/ and resources/. $QTPLANET_ROOT overrides the build-time default.
const std::string& assetRoot();

// Absolute path of an asset named relative to assetRoot().
std::string asset(const std::string& relative);

}   // namespace qtplanet

#endif   // QTPLANET_ASSETPATHS_H
