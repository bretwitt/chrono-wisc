#include "chrono_planet/lod/QuadtreeTile.h"

#include <cmath>

#include "chrono_planet/core/BenchProfiler.h"
#include "chrono_planet/dem/MultiGeoTIFFManager.h"
#include "chrono_planet/lod/SphericalCoordinates.h"

namespace {

// Deepest level the LOD builds.
constexpr int kMaxLeafLevel = 17;

// Children merge only once the camera has retreated this far past the split distance.
constexpr double kMergeHysteresis = 1.5;

// A half-built split is forgotten once the camera is this far past the merge distance.
constexpr double kForgetPendingFactor = 1.3;

// Tiles reaching this latitude are never refined.
constexpr double kPolarCapLatDeg = 89.9;

// True for a tile that touches the excluded polar cap.
template <typename Boundary>
bool inPolarCap(const Boundary& boundary) {
    return std::abs(boundary.centerLatDeg) + boundary.halfHeightDeg >= kPolarCapLatDeg;
}

}   // namespace

template <typename CoordSystem>
QuadtreeTile<CoordSystem>::QuadtreeTile(Boundary b, std::shared_ptr<const MultiGeoTIFFManager> geo)
    : geo_(std::move(geo)), tree_(std::make_unique<Node>(b)) {
    tree_->onInit = [this](Node* node) { buildMesh(node); };
    tree_->onDestroy = [this](Node* node) {
        pending_.erase(node);
        if (meshes_.erase(node)) {
            meshSetVersion_++;
        }
    };
    tree_->onSplit = [this](Node* parent) {
        parent->getType()->morphFromParent = true;
        for (Node* c : parent->children()) {
            c->getType()->morphFromParent = true;
        }
        meshSetVersion_++;
    };
    tree_->onMerge = [this](Node*) { meshSetVersion_++; };
    buildMesh(tree_.get());
}

template <typename CoordSystem>
QuadtreeTile<CoordSystem>::~QuadtreeTile() = default;

template <typename CoordSystem>
void QuadtreeTile<CoordSystem>::updateLOD(const qtplanet::Vec3& cameraM) {
    updateLODRec(tree_.get(), cameraM);
}

template <typename CoordSystem>
int QuadtreeTile<CoordSystem>::prefetch(const qtplanet::Vec3& predictedCameraM, Worker& worker, int budget) {
    return prefetchRec(tree_.get(), predictedCameraM, worker, budget);
}

template <typename CoordSystem>
void QuadtreeTile<CoordSystem>::takeBuilt(typename Worker::Result&& r) {
    // pending_ keys are live nodes (onDestroy erases them). The id guards against address reuse.
    auto it = pending_.find(static_cast<const Node*>(r.node));
    if (it == pending_.end() || it->first->getType()->nodeId != r.nodeId) {
        return;
    }
    PendingSplit& p = it->second;
    if (p.built[r.child]) {
        return;   // demand got there first
    }
    p.ready[r.child] = std::move(r.mesh);
    p.built[r.child] = true;
}

template <typename CoordSystem>
void QuadtreeTile<CoordSystem>::updateLODRec(Node* node, const qtplanet::Vec3& cameraM) {
    const Boundary boundary = node->getBoundary();
    if (inPolarCap(boundary)) {
        return;
    }
    const TileMetadata& meta = *node->getType();
    const double distance = CoordinateTraits<CoordSystem>::distanceToBounds(boundary, cameraM, meta.minElevation, meta.maxElevation);
    const int level = node->getLevel();
    const double splitDist = lodSplitDistance(level);
    const double mergeDist = splitDist * kMergeHysteresis;
    const bool canSplit = level < kMaxLeafLevel;

    if (canSplit && distance < splitDist) {
        if (!node->isDivided()) {
            // The split test decides. Whatever prefetch did not finish is built here, in this tick.
            BENCH_SCOPE("split");
            PendingSplit& p = completePending(node);
            inject_ = &p.ready;
            injectNext_ = 0;
            node->subdivide();   // children take their meshes from inject_
            inject_ = nullptr;
            pending_.erase(node);
            bench::bump(bench::ctr().splitsDone);
        }
        for (Node* c : node->children()) {
            updateLODRec(c, cameraM);
        }
    } else if (distance > mergeDist) {
        if (node->isDivided()) {
            BENCH_SCOPE("merge");
            node->getType()->morphFromParent = false;   // back at full detail, no split ease
            node->merge();
            bench::bump(bench::ctr().merges);
        }
        if (distance > mergeDist * kForgetPendingFactor) {
            pending_.erase(node);
        }
    }
}

template <typename CoordSystem>
int QuadtreeTile<CoordSystem>::prefetchRec(Node* node, const qtplanet::Vec3& predictedCameraM, Worker& worker, int budget) {
    int queued = 0;
    if (budget <= 0) {
        return queued;
    }
    const Boundary boundary = node->getBoundary();
    if (inPolarCap(boundary) || node->getLevel() >= kMaxLeafLevel) {
        return queued;
    }
    if (node->isDivided()) {
        for (Node* c : node->children()) {
            queued += prefetchRec(c, predictedCameraM, worker, budget - queued);
        }
        return queued;
    }
    // Same test updateLODRec will make, at the camera's predicted position.
    const TileMetadata& meta = *node->getType();
    const double distance = CoordinateTraits<CoordSystem>::distanceToBounds(boundary, predictedCameraM, meta.minElevation, meta.maxElevation);
    if (distance < lodSplitDistance(node->getLevel())) {
        queued += requestPending(node, worker, budget);
    }
    return queued;
}

template <typename CoordSystem>
int QuadtreeTile<CoordSystem>::requestPending(Node* node, Worker& worker, int budget) {
    int queued = 0;
    PendingSplit& p = pending_[node];
    for (int k = 0; k < 4 && queued < budget; ++k) {
        if (p.built[k] || p.requested[k]) {
            continue;
        }
        if (worker.inFlight() >= kMaxPrefetchInFlight) {
            return queued;
        }
        const auto bounds = CoordinateTraits<CoordSystem>::getChildBounds(node->getBoundary());
        worker.enqueue({this, node, node->getType()->nodeId, k, bounds[k], node->getLevel() + 1});
        p.requested[k] = true;
        ++queued;
    }
    return queued;
}

template <typename CoordSystem>
void QuadtreeTile<CoordSystem>::appendPendingMeshes(std::vector<const Mesh*>& out) const {
    for (const auto& kv : pending_) {
        for (int k = 0; k < 4; ++k) {
            if (kv.second.built[k]) {
                out.push_back(&kv.second.ready[k]);
            }
        }
    }
}

template <typename CoordSystem>
typename QuadtreeTile<CoordSystem>::PendingSplit& QuadtreeTile<CoordSystem>::completePending(Node* node) {
    PendingSplit& p = pending_[node];
    if (!p.allBuilt()) {
        // A child the worker is still on is built here anyway. Its late result is dropped by takeBuilt.
        const auto bounds = CoordinateTraits<CoordSystem>::getChildBounds(node->getBoundary());
        for (int k = 0; k < 4; ++k) {
            if (!p.built[k]) {
                p.ready[k] = buildTimed(bounds[k], node->getLevel() + 1);
                p.built[k] = true;
            }
        }
    }
    return p;
}

template <typename CoordSystem>
void QuadtreeTile<CoordSystem>::buildMesh(Node* node) {
    Mesh m = inject_ ? std::move((*inject_)[injectNext_++]) : buildTimed(node->getBoundary(), node->getLevel());
    TileMetadata& meta = *node->getType();
    meta.maxElevation = m.maxElevation;
    meta.minElevation = m.minElevation;
    meta.bake = std::move(m.bake);
    meshes_[node] = std::move(m);
    meshSetVersion_++;
}

template <typename CoordSystem>
Mesh QuadtreeTile<CoordSystem>::buildTimed(const Boundary& bounds, int level) const {
    BENCH_SCOPE("build_sync");
    const auto t0 = bench::Clock::now();
    Mesh m = generateMesh(bounds, level);
    bench::bump(bench::ctr().syncBuildNs,
                std::chrono::duration_cast<std::chrono::nanoseconds>(bench::Clock::now() - t0).count());
    return m;
}

template <typename CoordSystem>
Mesh QuadtreeTile<CoordSystem>::generateMesh(const Boundary& bounds, int level) const {
    return TileMeshBuilder<CoordSystem>(geo_).build(bounds, level);
}

template <typename CoordSystem>
double QuadtreeTile<CoordSystem>::getElevation(Position pos, int zoomLevel) const {
    return CoordinateTraits<CoordSystem>::computeBaseElevation(pos, geo_.get(), zoomLevel);
}

template class QuadtreeTile<Spherical>;
