#ifndef QTPLANET_QUADTREE_WORLD_H
#define QTPLANET_QUADTREE_WORLD_H

#include "chrono_planet/ChApiPlanet.h"

#include <algorithm>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "chrono_planet/core/MathUtil.h"
#include "chrono_planet/core/SphereMath.h"
#include "chrono_planet/dem/MultiGeoTIFFManager.h"
#include "chrono_planet/lod/QuadtreeTile.h"
#include "chrono_planet/lod/SphericalCoordinates.h"
#include "chrono_planet/lod/TileBuildWorker.h"

// Streams a ring of root tiles around the camera, each refined by its own LOD tree.
class CH_PLANET_API QuadtreeWorld {
public:
    // A DEM and the quadtree zoom range it provides elevation for.
    struct DemSource {
        std::string path;   // absolute, or relative to qtplanet::assetRoot()
        int minZoom;
        int maxZoom;
    };

    // The stock elevation stack, global LOLA at coarse zooms and the Apollo-region DEM below.
    static std::vector<DemSource> defaultDems();

    // Consumers that sample their own elevation must pass the same stack. Empty selects defaultDems().
    // Throws std::runtime_error when a DEM cannot be opened.
    QuadtreeWorld(double tileSizeDeg, int viewRangeInTiles, const std::vector<DemSource>& dems = {});
    ~QuadtreeWorld();
    QuadtreeWorld(const QuadtreeWorld&) = delete;
    QuadtreeWorld& operator=(const QuadtreeWorld&) = delete;

    // Sim-time length of one LOD tick. The tree at a sim time depends only on the camera path.
    static constexpr double kTickS = 1.0 / 60.0;

    // How far ahead of the camera, in sim seconds, tiles are prefetched.
    static constexpr double kPrefetchHorizonS = 1.0;

    // Camera position at sim time t.
    using PoseFn = std::function<qtplanet::Vec3(double t)>;

    // Advances the LOD to simTimeS, sampling poseAt at each tick along the way.
    void update(const PoseFn& poseAt, double simTimeS);
    // Same, for callers with only the current position. Ticks interpolate from the previous call.
    void update(double cameraX, double cameraY, double cameraZ, double simTimeS);

    // Child builds the prefetch worker may take per tick. Affects frame time only, never the tree.
    void setPrefetchBudget(int perTick) { prefetchBudget_ = std::max(0, perTick); }

    // Last tick stepped. A frame draws the tree this tick ended with.
    long long tick() const { return tick_; }

    // Resident root tiles.
    int getTotalTiles() const { return static_cast<int>(tiles_.size()); }

    // Every leaf mesh. Pointers stay valid until the next update().
    std::unordered_map<const QuadTree<TileMetadata, Spherical>*, const Mesh*> getAllMeshes() const;

    // Meshes prefetched for splits still to come, in Mesh::id order. Pointers stay valid until the next update().
    std::vector<const Mesh*> pendingMeshes() const;

    // Surface height at (lon, lat) as a mesh at zoom carries it, or nullopt outside the resident roots.
    // Zoom 0 is the bare DEM. Deeper zooms add the procedural relief the leaf meshes have.
    [[nodiscard]] std::optional<double> getElevation(double lon, double lat, int zoom = 0) const;

    // Changes whenever any tile's mesh set or the root ring changes.
    unsigned long long meshSetVersion() const;

    std::shared_ptr<const MultiGeoTIFFManager> geoManager() const { return geoMgr_; }

    // Camera of the last tick, in degrees.
    qtplanet::LonLat getCameraPositionLonLat() const;
    // Camera of the last tick, in meters above the sphere.
    double getCameraPositionElevation() const;

private:
    // Runs one LOD tick at a camera position in meters.
    void step(const qtplanet::Vec3& cameraM);

    std::shared_ptr<MultiGeoTIFFManager> geoMgr_;
    double tileSizeDeg_;       // root tile side, degrees
    int viewRangeInTiles_;     // ring half-width in tiles at ground level
    qtplanet::Vec3 cameraM_;   // position of the last tick
    bool haveTick_ = false;
    long long tick_ = 0;
    double prevT_ = 0;   // previous update() call, for interpolation
    qtplanet::Vec3 prevPositionM_;
    bool havePrevTick_ = false;
    qtplanet::Vec3 previousTickM_;   // position of the tick before last, for extrapolation
    int prefetchBudget_ = QuadtreeTile<Spherical>::kDefaultPrefetchBudget;
    std::unordered_map<TileKey, std::unique_ptr<QuadtreeTile<Spherical>>> tiles_;
    unsigned long long rootTopologyVersion_ = 0;   // bumped when a root tile is added or removed
    // Declared after tiles_ so the thread joins before any tile it built for is released.
    std::unique_ptr<TileBuildWorker<Spherical>> worker_;
};

#endif   // QTPLANET_QUADTREE_WORLD_H
