#ifndef QTPLANET_QUADTREE_TILE_H
#define QTPLANET_QUADTREE_TILE_H

#include "chrono_planet/ChApiPlanet.h"

#include <algorithm>
#include <array>
#include <memory>
#include <unordered_map>
#include "chrono_planet/core/MathUtil.h"
#include "chrono_planet/core/Planet.h"
#include "chrono_planet/lod/Quadtree.h"
#include "chrono_planet/lod/TileBuildWorker.h"
#include "chrono_planet/lod/TileMeshBuilder.h"
#include "chrono_planet/lod/TileMetadata.h"

class MultiGeoTIFFManager;

// Camera distance below which a tile at level splits.
inline double lodSplitDistance(int level) {
    if (level < 0) {
        return qtplanet::kRadiusM;
    }
    return qtplanet::kRadiusM * 0.5 / static_cast<double>(1 << std::min(level, 30));
}

// One root tile with its LOD tree, each node's mesh, and child meshes prefetched for splits.
template <typename CoordSystem>
class QuadtreeTile {
public:
    using Boundary = typename CoordSystem::Boundary;
    using Position = typename CoordSystem::Position;
    using Node = QuadTree<TileMetadata, CoordSystem>;

    // Builds the root mesh synchronously from the DEM stack geo.
    QuadtreeTile(Boundary b, std::shared_ptr<const MultiGeoTIFFManager> geo);
    ~QuadtreeTile();
    QuadtreeTile(const QuadtreeTile&) = delete;
    QuadtreeTile& operator=(const QuadtreeTile&) = delete;

    // Child builds handed to the worker per tick. Affects speed only, never visibility.
    static constexpr int kDefaultPrefetchBudget = 2;
    // Jobs the worker may hold at once.
    static constexpr int kMaxPrefetchInFlight = 16;

    using Worker = TileBuildWorker<CoordSystem>;

    // Consumes a finished prefetch build for its node. Stale (merged or reused) nodes are dropped by id.
    void takeBuilt(typename Worker::Result&& r);
    // One LOD step at the tick's camera position (meters). Splits (building any missing child synchronously) and merges.
    void updateLOD(const qtplanet::Vec3& cameraM);
    // Queues up to budget child builds for leaves that will split at the predicted camera. Returns the number queued.
    [[nodiscard]] int prefetch(const qtplanet::Vec3& predictedCameraM, Worker& worker, int budget);

    // Every node's mesh, divided nodes included.
    const std::unordered_map<const Node*, Mesh>& getMeshes() const { return meshes_; }

    // Appends meshes built for splits not yet taken to out.
    void appendPendingMeshes(std::vector<const Mesh*>& out) const;

    // Changes whenever a mesh is added, regenerated or removed.
    unsigned meshSetVersion() const { return meshSetVersion_; }

    // Ground height at pos as a mesh at zoomLevel would carry it.
    [[nodiscard]] double getElevation(Position pos, int zoomLevel) const;

private:
    // Child meshes built ahead of a split, in getChildBounds order. Entry k of each flag set tracks child k.
    struct PendingSplit {
        std::array<Mesh, 4> ready;
        std::array<bool, 4> built{}, requested{};
        bool allBuilt() const {
            return std::all_of(built.begin(), built.end(), [](bool b) { return b; });
        }
    };

    void updateLODRec(Node* node, const qtplanet::Vec3& cameraM);
    int prefetchRec(Node* node, const qtplanet::Vec3& predictedCameraM, Worker& worker, int budget);
    int requestPending(Node* node, Worker& worker, int budget);   // queues one child build, within budget
    PendingSplit& completePending(Node* node);                    // builds every missing child now
    void buildMesh(Node* node);                                   // ensures the node has a mesh and metadata
    Mesh buildTimed(const Boundary& bounds, int level) const;     // generateMesh, timed for the bench
    Mesh generateMesh(const Boundary& bounds, int level) const;

    std::shared_ptr<const MultiGeoTIFFManager> geo_;
    std::unordered_map<const Node*, Mesh> meshes_;
    std::unordered_map<const Node*, PendingSplit> pending_;   // half-built splits, keyed by the leaf
    // During subdivide(), the meshes the new children take, in creation order. Borrowed from pending_.
    std::array<Mesh, 4>* inject_ = nullptr;
    int injectNext_ = 0;
    unsigned meshSetVersion_ = 0;
    // The root node. Declared last so its onDestroy callbacks still find the maps above alive.
    std::unique_ptr<Node> tree_;
};

#endif   // QTPLANET_QUADTREE_TILE_H
