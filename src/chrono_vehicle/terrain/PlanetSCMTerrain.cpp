// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2026 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: bgwitt
// =============================================================================

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "chrono_vehicle/terrain/PlanetSCMTerrain.h"

namespace chrono {
namespace vehicle {

// -----------------------------------------------------------------------------
// PlanetSCMHeightFunctor
// -----------------------------------------------------------------------------

PlanetSCMHeightFunctor::PlanetSCMHeightFunctor(std::shared_ptr<const planet::ChPlanetSurface> surface,
                                               const planet::ChSiteFrame& site)
    : m_surface(std::move(surface)), m_site(site) {}

double PlanetSCMHeightFunctor::GetInitHeight(const ChVector2i& loc, double delta) {
    return HeightOf(loc, delta);
}

double PlanetSCMHeightFunctor::HeightOf(const ChVector2i& loc, double delta) const {
    const std::int64_t k = Key(loc);
    Shard& shard = m_shards[ShardOf(k)];
    {
        std::lock_guard<std::mutex> lock(shard.mutex);
        auto it = shard.heights.find(k);
        if (it != shard.heights.end())
            return it->second;
    }
    // Evaluate outside the lock: the surface sample is the expensive part and is thread-safe.
    double lon, lat;
    m_site.ToLonLat(loc.x() * delta, loc.y() * delta, lon, lat);
    const double h = m_surface->GetElevation(lon, lat) - m_site.GetOriginElevation();
    std::lock_guard<std::mutex> lock(shard.mutex);
    shard.heights.emplace(k, h);
    return h;
}

void PlanetSCMHeightFunctor::GetInitHeights(const std::vector<ChVector2i>& locs,
                                            double delta,
                                            std::vector<double>& out) const {
    out.resize(locs.size());
    for (std::size_t i = 0; i < locs.size(); ++i)
        out[i] = HeightOf(locs[i], delta);
}

std::size_t PlanetSCMHeightFunctor::GetNumCachedNodes() const {
    std::size_t n = 0;
    for (const auto& shard : m_shards) {
        std::lock_guard<std::mutex> lock(shard.mutex);
        n += shard.heights.size();
    }
    return n;
}

// -----------------------------------------------------------------------------
// PlanetSCMTerrain
// -----------------------------------------------------------------------------

PlanetSCMTerrain::PlanetSCMTerrain(ChSystem* system,
                                   std::shared_ptr<const planet::ChPlanetSurface> surface,
                                   const planet::ChSiteFrame& site)
    : SCMTerrain(system, false), m_functor(std::make_shared<PlanetSCMHeightFunctor>(std::move(surface), site)) {}

PlanetSCMTerrain::~PlanetSCMTerrain() {
    m_stop.store(true, std::memory_order_release);
    // Notify under each mutex so a worker between its predicate test and wait cannot miss the wake.
    {
        std::lock_guard<std::mutex> lock(m_stats_mutex);
        m_stats_cv.notify_all();
    }
    {
        std::lock_guard<std::mutex> lock(m_prefetch_mutex);
        m_prefetch_cv.notify_all();
    }
    if (m_stats_thread.joinable())
        m_stats_thread.join();
    if (m_prefetch_thread.joinable())
        m_prefetch_thread.join();
}

void PlanetSCMTerrain::Initialize(const Params& params, const std::vector<Wheel>& wheels) {
    if (wheels.empty())
        throw std::invalid_argument(
            "PlanetSCMTerrain::Initialize: no wheels; SCM would ray-cast the AABB of every collision shape.");

    m_params = params;

    SetSoilParameters(m_params.bekker_kphi, m_params.bekker_kc, m_params.bekker_n, m_params.mohr_cohesion,
                      m_params.mohr_friction, m_params.janosi_shear, m_params.elastic_k, m_params.damping_r);
    SetTestHeight(m_params.test_height);
    if (m_params.bulldozing) {
        EnableBulldozing(true);
        SetBulldozingParameters(m_params.erosion_angle);
    }

    // Active domains are registered before Initialize. The box is in the wheel's reference frame with full
    // dimensions: the wheel spins about its local y, so x and z span the diameter and y the width.
    for (const auto& w : wheels) {
        const double d = 2.0 * w.radius + 2.0 * m_params.domain_pad;
        const double b = w.width + 2.0 * m_params.domain_pad;
        AddActiveDomain(w.body, ChVector3d(0, 0, 0), ChVector3d(d, b, d));
    }

    SCMTerrain::Initialize(m_functor, m_params.delta);

    // Workers start only once the grid exists.
    m_stats_thread = std::thread([this] { StatsLoop(); });
    if (m_params.prefetch)
        m_prefetch_thread = std::thread([this] { PrefetchLoop(); });
}

void PlanetSCMTerrain::SetSoil(const Params& params) {
    // Grid geometry and the prefetch worker are fixed at Initialize; keep ours.
    const double delta = m_params.delta;
    const double domain_pad = m_params.domain_pad;
    const bool prefetch = m_params.prefetch;
    const double lookahead = m_params.prefetch_lookahead;
    const double half_width = m_params.prefetch_half_width;
    m_params = params;
    m_params.delta = delta;
    m_params.domain_pad = domain_pad;
    m_params.prefetch = prefetch;
    m_params.prefetch_lookahead = lookahead;
    m_params.prefetch_half_width = half_width;

    SetSoilParameters(m_params.bekker_kphi, m_params.bekker_kc, m_params.bekker_n, m_params.mohr_cohesion,
                      m_params.mohr_friction, m_params.janosi_shear, m_params.elastic_k, m_params.damping_r);
    SetTestHeight(m_params.test_height);
    EnableBulldozing(m_params.bulldozing);
    if (m_params.bulldozing)
        SetBulldozingParameters(m_params.erosion_angle);
}

double PlanetSCMTerrain::GetContactForce(const std::shared_ptr<ChBody>& body) const {
    ChVector3d force, torque;
    if (!GetContactForceBody(body, force, torque))
        return 0.0;
    return force.Length();
}

// -----------------------------------------------------------------------------
// Rut statistics. The walk grows with the ground driven over, so only the snapshot happens on the
// simulation thread (Chrono's node map is read between steps); the arithmetic runs on the worker.
// -----------------------------------------------------------------------------

void PlanetSCMTerrain::RequestRutStats() {
    {
        std::lock_guard<std::mutex> lock(m_stats_mutex);
        if (m_stats_in_flight)
            return;
        m_stats_in_flight = true;
    }

    auto nodes = GetModifiedNodes(true);

    std::lock_guard<std::mutex> lock(m_stats_mutex);
    if (nodes.empty()) {
        // An empty job would never wake the worker and would leave the in-flight flag stuck.
        m_stats_result = RutStats{};
        m_stats_ready = true;
        m_stats_in_flight = false;
        return;
    }
    m_stats_job = std::move(nodes);
    m_stats_cv.notify_one();
}

std::optional<PlanetSCMTerrain::RutStats> PlanetSCMTerrain::TakeRutStats() {
    std::lock_guard<std::mutex> lock(m_stats_mutex);
    if (!m_stats_ready)
        return std::nullopt;
    m_stats_ready = false;
    return m_stats_result;
}

void PlanetSCMTerrain::StatsLoop() {
    std::vector<SCMTerrain::NodeLevel> job;
    std::vector<ChVector2i> locs;
    std::vector<double> init;
    const double delta = m_params.delta;

    while (true) {
        {
            std::unique_lock<std::mutex> lock(m_stats_mutex);
            m_stats_cv.wait(lock, [this] { return !m_stats_job.empty() || m_stop.load(std::memory_order_acquire); });
            if (m_stop.load(std::memory_order_acquire))
                return;
            job.swap(m_stats_job);
            m_stats_job.clear();
        }

        RutStats st;
        st.nodes = job.size();
        st.area_m2 = static_cast<double>(st.nodes) * delta * delta;

        locs.clear();
        locs.reserve(job.size());
        for (const auto& nl : job)
            locs.push_back(nl.first);
        m_functor->GetInitHeights(locs, delta, init);

        double sum = 0.0;
        for (std::size_t v = 0; v < job.size(); ++v) {
            const double d = init[v] - job[v].second;  // undeformed minus current level
            if (d <= 0.0)
                continue;  // untouched ground, or a bulldozed rim
            sum += d;
            ++st.deformed;
            if (d > st.max_depth_m)
                st.max_depth_m = d;
        }
        st.mean_depth_m = st.deformed ? sum / st.deformed : 0.0;
        job.clear();

        std::lock_guard<std::mutex> lock(m_stats_mutex);
        m_stats_result = st;
        m_stats_ready = true;
        m_stats_in_flight = false;
    }
}

// -----------------------------------------------------------------------------
// Undeformed-height prefetch
// -----------------------------------------------------------------------------

namespace {
// Floor division, so a corridor crossing the site origin does not fold negative blocks onto positive ones.
int FloorDiv(int a, int b) {
    const int q = a / b;
    return (a % b != 0 && ((a < 0) != (b < 0))) ? q - 1 : q;
}
}  // namespace

void PlanetSCMTerrain::SetPrefetchPose(double x, double y, double hx, double hy) {
    if (!m_params.prefetch)
        return;

    const double n = std::hypot(hx, hy);
    if (n < 1e-9) {
        hx = 1.0;
        hy = 0.0;
    } else {
        hx /= n;
        hy /= n;
    }

    std::lock_guard<std::mutex> lock(m_prefetch_mutex);
    // Only wake the worker once the corridor has actually moved (a quarter meter or about 8 degrees).
    const bool moved = std::hypot(x - m_pre_x, y - m_pre_y) > 0.25 || (hx * m_pre_hx + hy * m_pre_hy) < 0.99;
    if (m_pre_notified && !moved)
        return;
    m_pre_x = x;
    m_pre_y = y;
    m_pre_hx = hx;
    m_pre_hy = hy;
    m_pre_dirty = true;
    m_pre_notified = true;
    m_prefetch_cv.notify_one();
}

void PlanetSCMTerrain::PrefetchLoop() {
    while (true) {
        double x, y, hx, hy;
        {
            std::unique_lock<std::mutex> lock(m_prefetch_mutex);
            m_prefetch_cv.wait(lock, [this] { return m_pre_dirty || m_stop.load(std::memory_order_acquire); });
            if (m_stop.load(std::memory_order_acquire))
                return;
            m_pre_dirty = false;
            x = m_pre_x;
            y = m_pre_y;
            hx = m_pre_hx;
            hy = m_pre_hy;
        }

        if (!WarmCorridor(x, y, hx, hy)) {
            // Stopped on the per-pass cap with blocks left: re-arm rather than wait for the next pose.
            std::lock_guard<std::mutex> lock(m_prefetch_mutex);
            m_pre_dirty = true;
        }
    }
}

bool PlanetSCMTerrain::WarmCorridor(double x, double y, double hx, double hy) {
    const double delta = m_params.delta;
    const double lookahead = m_params.prefetch_lookahead;
    const double half_width = m_params.prefetch_half_width;

    // A block is a square meter of ground, so the warmed set stays small over a long drive.
    const int nb = std::max(1, static_cast<int>(std::lround(1.0 / delta)));
    const double block_m = nb * delta;
    const double pad = block_m * 0.70710678;  // block half-diagonal

    // Half a width behind as well as ahead, so turning on the spot still lands on warm ground.
    const double back = half_width;
    double min_x = x, max_x = x, min_y = y, max_y = y;
    for (const double a : {-back, lookahead}) {
        for (const double l : {-half_width, half_width}) {
            const double px = x + a * hx - l * hy;
            const double py = y + a * hy + l * hx;
            min_x = std::min(min_x, px);
            max_x = std::max(max_x, px);
            min_y = std::min(min_y, py);
            max_y = std::max(max_y, py);
        }
    }

    const int bi0 = FloorDiv(static_cast<int>(std::floor((min_x - pad) / delta)), nb);
    const int bi1 = FloorDiv(static_cast<int>(std::ceil((max_x + pad) / delta)), nb);
    const int bj0 = FloorDiv(static_cast<int>(std::floor((min_y - pad) / delta)), nb);
    const int bj1 = FloorDiv(static_cast<int>(std::ceil((max_y + pad) / delta)), nb);

    // Bound the pass so a pose jump cannot hold shutdown on the join.
    constexpr int kMaxBlocksPerPass = 512;
    int warmed = 0;

    for (int bj = bj0; bj <= bj1; ++bj) {
        for (int bi = bi0; bi <= bi1; ++bi) {
            if (m_stop.load(std::memory_order_acquire))
                return true;
            if (warmed >= kMaxBlocksPerPass)
                return false;

            const std::int64_t k = (static_cast<std::int64_t>(bi) << 32) | static_cast<std::uint32_t>(bj);
            if (m_warmed_blocks.count(k))
                continue;

            // Test the block center against the corridor grown by the block half-diagonal.
            const double cx = (bi * nb + nb * 0.5) * delta;
            const double cy = (bj * nb + nb * 0.5) * delta;
            const double dx = cx - x, dy = cy - y;
            const double along = dx * hx + dy * hy;
            const double lat = -dx * hy + dy * hx;
            if (along < -back - pad || along > lookahead + pad)
                continue;
            if (std::fabs(lat) > half_width + pad)
                continue;

            for (int j = bj * nb; j < (bj + 1) * nb; ++j)
                for (int i = bi * nb; i < (bi + 1) * nb; ++i)
                    m_functor->Warm(ChVector2i(i, j), delta);
            m_warmed_blocks.insert(k);
            ++warmed;
        }
    }
    return true;
}

}  // end namespace vehicle
}  // end namespace chrono
