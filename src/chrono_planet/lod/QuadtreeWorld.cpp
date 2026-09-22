#include "chrono_planet/lod/QuadtreeWorld.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_set>

#include "chrono_planet/core/AssetPaths.h"

using Traits = CoordinateTraits<Spherical>;

namespace {

// The active ring widens by one tile per this much altitude.
constexpr double kAltitudePerExtraTileM = 10000.0;
// Ring half-width bounds in tiles.
constexpr int kMinRingTiles = 2, kMaxRingTiles = 90;
// Keeps root add/remove events distinct from per-tile mesh changes in meshSetVersion().
constexpr unsigned long long kRootVersionStride = 1000003ULL;

bool isAbsolute(const std::string& path) {
    return !path.empty() && path.front() == '/';
}

}   // namespace

std::vector<QuadtreeWorld::DemSource> QuadtreeWorld::defaultDems() {
    return {
        {"resources/ldem_64_fixed.tif", 0, 4},
        {"resources/ldem_1024_apollo_region.tif", 5, 30},
    };
}

QuadtreeWorld::QuadtreeWorld(double tileSizeDeg, int viewRangeInTiles, const std::vector<DemSource>& dems)
    : geoMgr_(std::make_shared<MultiGeoTIFFManager>()), tileSizeDeg_(tileSizeDeg), viewRangeInTiles_(viewRangeInTiles) {
    const std::vector<DemSource> sources = dems.empty() ? defaultDems() : dems;
    for (const DemSource& d : sources) {
        const std::string path = isAbsolute(d.path) ? d.path : qtplanet::asset(d.path);
        if (!geoMgr_->addSource(path, d.minZoom, d.maxZoom)) {
            throw std::runtime_error("QuadtreeWorld: cannot open DEM " + path);
        }
    }
    worker_ = std::make_unique<TileBuildWorker<Spherical>>(geoMgr_);
}

QuadtreeWorld::~QuadtreeWorld() = default;

void QuadtreeWorld::update(const PoseFn& poseAt, double simTimeS) {
    // The tick this time falls in. The epsilon keeps t = k / 60 exactly on tick k.
    const long long tick = static_cast<long long>(std::floor(simTimeS / kTickS + 1e-6));
    const bool restart = !haveTick_ || simTimeS < prevT_;   // first call, or the sim clock was reset
    if (restart) {
        haveTick_ = true;
        havePrevTick_ = false;
    }
    // Step every tick since the last call, or only this one after a restart.
    for (long long k = restart ? tick : tick_ + 1; k <= tick; ++k) {
        step(poseAt(k * kTickS));
    }
    tick_ = tick;
    prevT_ = simTimeS;
}

void QuadtreeWorld::update(double camX, double camY, double camZ, double simTimeS) {
    const double t0 = prevT_;
    const qtplanet::Vec3 from = prevPositionM_, to{camX, camY, camZ};
    const bool have = haveTick_ && simTimeS > t0;
    update(
        [&](double t) {
            const double u = have ? std::clamp((t - t0) / (simTimeS - t0), 0.0, 1.0) : 1.0;
            return qtplanet::Vec3{qtplanet::lerp(from.x, to.x, u), qtplanet::lerp(from.y, to.y, u), qtplanet::lerp(from.z, to.z, u)};
        },
        simTimeS);
    prevPositionM_ = to;
}

void QuadtreeWorld::step(const qtplanet::Vec3& cameraM) {
    // The camera kPrefetchHorizonS ahead.
    qtplanet::Vec3 predictedM = cameraM;
    if (havePrevTick_) {
        const double ahead = kPrefetchHorizonS / kTickS;
        predictedM.x += (cameraM.x - previousTickM_.x) * ahead;
        predictedM.y += (cameraM.y - previousTickM_.y) * ahead;
        predictedM.z += (cameraM.z - previousTickM_.z) * ahead;
    }
    havePrevTick_ = true;
    previousTickM_ = cameraM_ = cameraM;

    const qtplanet::LonLat cam = qtplanet::lonLatOf(cameraM.x, cameraM.y, cameraM.z);
    const auto [centreX, centreY] = Traits::computeTileIndices({cam.lon, cam.lat}, tileSizeDeg_);

    // Widen the active ring with altitude.
    const int extraTiles = static_cast<int>(std::ceil(Traits::elevationOf(cameraM.x, cameraM.y, cameraM.z) / kAltitudePerExtraTileM));
    const int range = std::clamp(viewRangeInTiles_ + extraTiles, kMinRingTiles, kMaxRingTiles);

    const int nLon = static_cast<int>(std::round(360.0 / tileSizeDeg_));
    const int nLat = static_cast<int>(std::round(180.0 / tileSizeDeg_));

    std::unordered_set<TileKey> needed;
    for (int dy = -range; dy <= range; ++dy) {
        for (int dx = -range; dx <= range; ++dx) {
            const TileKey key{((centreX + dx) % nLon + nLon) % nLon, std::clamp(centreY + dy, 0, nLat - 1)};
            needed.insert(key);
            if (tiles_.count(key)) {
                continue;
            }
            const auto [cLon, cLat] = Traits::tileCenterPosition(key, tileSizeDeg_);
            const double half = tileSizeDeg_ * 0.5;
            tiles_[key] = std::make_unique<QuadtreeTile<Spherical>>(Spherical::Boundary{cLon, cLat, half, half}, geoMgr_);
            rootTopologyVersion_++;
        }
    }

    for (auto it = tiles_.begin(); it != tiles_.end();) {
        if (needed.count(it->first)) {
            ++it;
            continue;
        }
        it = tiles_.erase(it);
        rootTopologyVersion_++;
    }

    // Hand the worker's finished builds to the tiles that asked. A root that left the ring is dropped.
    for (auto& r : worker_->drain()) {
        for (auto& p : tiles_) {
            if (p.second.get() == r.owner) {
                p.second->takeBuilt(std::move(r));
                break;
            }
        }
    }
    for (auto& p : tiles_) {
        p.second->updateLOD(cameraM);
    }
    int queued = 0;            // count of requests made this tick
    for (auto& p : tiles_) {   // children the camera is about to need
        queued += p.second->prefetch(predictedM, *worker_, prefetchBudget_ - queued);
    }
}

std::unordered_map<const QuadTree<TileMetadata, Spherical>*, const Mesh*> QuadtreeWorld::getAllMeshes() const {
    std::unordered_map<const QuadTree<TileMetadata, Spherical>*, const Mesh*> out;
    for (const auto& p : tiles_) {
        for (const auto& kv : p.second->getMeshes()) {
            if (!kv.first->isDivided()) {
                out.emplace(kv.first, &kv.second);   // divided nodes keep a hidden mesh
            }
        }
    }
    return out;
}

std::vector<const Mesh*> QuadtreeWorld::pendingMeshes() const {
    std::vector<const Mesh*> out;
    for (const auto& p : tiles_) {
        p.second->appendPendingMeshes(out);
    }
    // Build order, so a per-frame upload cap reaches the same meshes first in every run.
    std::sort(out.begin(), out.end(), [](const Mesh* a, const Mesh* b) { return a->id < b->id; });
    return out;
}

std::optional<double> QuadtreeWorld::getElevation(double lon, double lat, int zoom) const {
    const auto [x, y] = Traits::computeTileIndices({lon, lat}, tileSizeDeg_);
    const auto it = tiles_.find(TileKey{x, y});
    if (it == tiles_.end()) {
        return std::nullopt;
    }
    return it->second->getElevation({lon, lat}, zoom);
}

unsigned long long QuadtreeWorld::meshSetVersion() const {
    unsigned long long v = rootTopologyVersion_ * kRootVersionStride;
    for (const auto& p : tiles_) {
        v += p.second->meshSetVersion();
    }
    return v;
}

qtplanet::LonLat QuadtreeWorld::getCameraPositionLonLat() const {
    return qtplanet::lonLatOf(cameraM_.x, cameraM_.y, cameraM_.z);
}

double QuadtreeWorld::getCameraPositionElevation() const {
    return Traits::elevationOf(cameraM_.x, cameraM_.y, cameraM_.z);
}
