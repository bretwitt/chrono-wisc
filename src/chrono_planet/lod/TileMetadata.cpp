#include "chrono_planet/lod/TileMetadata.h"

#include <atomic>

std::uint64_t ShadingBake::nextId() {
    static std::atomic<std::uint64_t> n{1};
    return n++;
}

std::uint64_t TileMetadata::nextNodeId() {
    static std::atomic<std::uint64_t> n{1};
    return n++;
}
