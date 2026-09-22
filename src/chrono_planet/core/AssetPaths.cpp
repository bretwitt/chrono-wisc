#include "chrono_planet/core/AssetPaths.h"

#include <cstdlib>
#include <sys/stat.h>

namespace qtplanet {
namespace {

bool hasAssets(const std::string& dir) {
    struct stat st;
    return !dir.empty() && stat((dir + "/shaders").c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

// Three levels up from this file, which the source list records by absolute path.
std::string compiledFromRoot() {
    const std::string self = __FILE__;
    if (self.empty() || self[0] != '/') {
        return {};
    }
    auto cut = [&](std::string s) {
        const auto slash = s.find_last_of('/');
        return slash == std::string::npos ? std::string{} : s.substr(0, slash);
    };
    return cut(cut(cut(cut(self))));   // .../src/qtplanet/core/AssetPaths.cpp -> ...
}

std::string resolve() {
    if (const char* env = std::getenv("QTPLANET_ROOT"); env && *env) {
        return env;
    }
#ifdef QTPLANET_SOURCE_ROOT
    if (hasAssets(QTPLANET_SOURCE_ROOT)) {
        return QTPLANET_SOURCE_ROOT;
    }
#endif
    if (const std::string root = compiledFromRoot(); hasAssets(root)) {
        return root;
    }
    return "..";
}

}   // namespace

const std::string& assetRoot() {
    static const std::string root = resolve();
    return root;
}

std::string asset(const std::string& relative) {
    return assetRoot() + "/" + relative;
}

}   // namespace qtplanet
