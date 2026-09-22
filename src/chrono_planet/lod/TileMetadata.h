#ifndef QTPLANET_TILEMETADATA_H
#define QTPLANET_TILEMETADATA_H

#include "chrono_planet/ChApiPlanet.h"

#include <cstdint>
#include <memory>
#include <vector>

// Identifies a root tile by its grid indices.
struct CH_PLANET_API TileKey {
    int x;
    int y;
    bool operator==(const TileKey& other) const { return x == other.x && y == other.y; }
};

namespace std {
// Hash for TileKey, so it can key QuadtreeWorld's root tile map.
template <>
struct hash<TileKey> {
    std::size_t operator()(const TileKey& key) const {
        return std::hash<int>()(key.x) ^ (std::hash<int>()(key.y) << 1);
    }
};
}   // namespace std

// Per-tile normal and height maps baked from the DEM stack, keyed on the GPU by their process-unique id.
struct CH_PLANET_API ShadingBake {
    static std::uint64_t nextId();
    std::uint64_t id = nextId();
    // The normal map covers the tile plus Mesh::kBakeMargin.
    int res = 0;                                             // texels per side
    double minLon = 0, maxLon = 0, minLat = 0, maxLat = 0;   // covered rect (deg), incl. margin
    std::vector<unsigned short> normals;                     // snorm16 xyz per texel, row-major, row 0 = minLat
    // The height map covers the tile plus Mesh::kHeightMargin, on its own coarser grid.
    int hRes = 0;   // texels per side
    double hMinLon = 0, hMaxLon = 0, hMinLat = 0, hMaxLat = 0;
    std::vector<float> heights;   // row-major, row 0 = hMinLat
};

// Per-node LOD state, stored as the payload of each QuadTree node.
struct CH_PLANET_API TileMetadata {
    // Process-unique per node, so a build for a node that has since merged is dropped.
    static std::uint64_t nextNodeId();
    std::uint64_t nodeId = nextNodeId();
    // Elevation range of the node's mesh (m above the sphere), for the LOD distance test.
    double maxElevation = 0.0;
    double minElevation = 0.0;
    bool morphFromParent = false;              // new leaves start their renderer-owned split animation
    std::shared_ptr<const ShadingBake> bake;   // null below Mesh::kBakeMinLevel
};

#endif   // QTPLANET_TILEMETADATA_H
