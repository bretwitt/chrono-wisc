#include "chrono_planet/lod/TileBuildWorker.h"

#include <chrono>

#include "chrono_planet/core/BenchProfiler.h"
#include "chrono_planet/lod/SphericalCoordinates.h"

template <typename C>
TileBuildWorker<C>::TileBuildWorker(std::shared_ptr<const MultiGeoTIFFManager> geo)
    : geo_(std::move(geo)), thread_([this]() { run(); }) {}

template <typename C>
TileBuildWorker<C>::~TileBuildWorker() {
    {
        std::lock_guard<std::mutex> lock(m_);
        stop_ = true;
        jobs_.clear();
    }
    cv_.notify_all();
    if (thread_.joinable()) {
        thread_.join();
    }
}

template <typename C>
void TileBuildWorker<C>::enqueue(Job job) {
    {
        std::lock_guard<std::mutex> lock(m_);
        jobs_.push_back(std::move(job));
    }
    cv_.notify_one();
}

template <typename C>
int TileBuildWorker<C>::inFlight() const {
    std::lock_guard<std::mutex> lock(m_);
    return static_cast<int>(jobs_.size()) + building_;
}

template <typename C>
std::vector<typename TileBuildWorker<C>::Result> TileBuildWorker<C>::drain() {
    std::vector<Result> out;
    std::lock_guard<std::mutex> lock(m_);
    out.swap(results_);
    return out;
}

template <typename C>
void TileBuildWorker<C>::run() {
    bench::Prof::get().setLabel("worker");
    for (;;) {
        Job job;
        {
            std::unique_lock<std::mutex> lock(m_);
            cv_.wait(lock, [this]() { return stop_ || !jobs_.empty(); });
            if (stop_) {
                return;
            }
            job = std::move(jobs_.front());
            jobs_.pop_front();
            ++building_;
        }
        bench::bump(bench::ctr().jobsStarted);
        bench::bump(bench::ctr().jobsInFlight);
        const auto t0 = bench::Clock::now();
        Mesh mesh = TileMeshBuilder<C>(geo_).build(job.bounds, job.level);
        const long long ns = std::chrono::duration_cast<std::chrono::nanoseconds>(bench::Clock::now() - t0).count();
        bench::bump(bench::ctr().jobNsTotal, ns);
        bench::bumpMax(bench::ctr().jobNsMax, ns);
        bench::bump(bench::ctr().jobsFinished);
        bench::bump(bench::ctr().jobsInFlight, -1);
        {
            std::lock_guard<std::mutex> lock(m_);
            --building_;
            results_.push_back(Result{job.owner, job.node, job.nodeId, job.child, std::move(mesh)});
        }
    }
}

template class TileBuildWorker<Spherical>;
