#ifndef QTPLANET_TILEBUILDWORKER_H
#define QTPLANET_TILEBUILDWORKER_H

#include "chrono_planet/ChApiPlanet.h"

#include <condition_variable>
#include <cstdint>
#include <deque>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#include "chrono_planet/lod/TileMeshBuilder.h"

class MultiGeoTIFFManager;

// Builds prefetched child meshes on its own thread. Readiness never gates visibility.
template <typename CoordSystem>
class TileBuildWorker {
public:
    using Boundary = typename CoordSystem::Boundary;

    // One child mesh to build. Its keys are echoed in the Result so the owner can route it.
    struct Job {
        const void* owner;      // the QuadtreeTile, as an opaque key
        const void* node;       // the node the children belong to, likewise
        std::uint64_t nodeId;   // TileMetadata::nodeId, against address reuse
        int child;              // 0..3, Traits::getChildBounds order
        Boundary bounds;
        int level;
    };
    // A finished Job with its mesh.
    struct Result {
        const void* owner;
        const void* node;
        std::uint64_t nodeId;
        int child;
        Mesh mesh;
    };

    // Starts the thread.
    explicit TileBuildWorker(std::shared_ptr<const MultiGeoTIFFManager> geo);
    // Drops queued jobs, finishes the one in progress and joins.
    ~TileBuildWorker();
    TileBuildWorker(const TileBuildWorker&) = delete;
    TileBuildWorker& operator=(const TileBuildWorker&) = delete;

    // Queues a build. Never blocks.
    void enqueue(Job job);
    // Jobs queued or being built.
    int inFlight() const;
    // Every finished build since the last call.
    [[nodiscard]] std::vector<Result> drain();

private:
    // Thread body.
    void run();

    std::shared_ptr<const MultiGeoTIFFManager> geo_;
    mutable std::mutex m_;
    std::condition_variable cv_;
    std::deque<Job> jobs_;
    std::vector<Result> results_;
    int building_ = 0;   // jobs taken off the queue and not yet in results_
    bool stop_ = false;
    std::thread thread_;
};

#endif   // QTPLANET_TILEBUILDWORKER_H
