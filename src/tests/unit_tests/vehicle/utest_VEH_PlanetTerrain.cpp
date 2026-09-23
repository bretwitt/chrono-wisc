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
//
// Unit tests for the planetary terrain models. No DEM files: the surface is the
// procedural fallback plus craters, which is enough to check the frame mapping,
// the height and normal queries, patch following and the SCM height functor.
//
// =============================================================================

#include <cmath>
#include <chrono>
#include <optional>
#include <thread>

#include "gtest/gtest.h"

#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemNSC.h"

#include "chrono_planet/ChPlanetSurface.h"
#include "chrono_planet/ChSiteFrame.h"
#include "chrono_planet/planets/moon/ChMoon.h"

#include "chrono_vehicle/terrain/PlanetSCMTerrain.h"
#include "chrono_vehicle/terrain/PlanetTerrain.h"

using namespace chrono;
using namespace chrono::planet;
using namespace chrono::vehicle;

namespace {

const double site_lon = 30.75;
const double site_lat = 20.19;
const int zoom = 15;

std::shared_ptr<ChPlanetSurface> MakeSurface() {
    return moon::CreateSurface(zoom);
}

}  // namespace

TEST(ChSiteFrame, RoundTrip) {
    ChSiteFrame site(moon::Body(), site_lon, site_lat, 100.0);
    for (double x : {-500.0, 0.0, 250.0}) {
        for (double y : {-300.0, 0.0, 700.0}) {
            double lon, lat;
            site.ToLonLat(x, y, lon, lat);
            const ChVector3d p = site.ToLocal(lon, lat, 100.0 + 5.0);
            EXPECT_NEAR(p.x(), x, 1e-6);
            EXPECT_NEAR(p.y(), y, 1e-6);
            EXPECT_NEAR(p.z(), 5.0, 1e-9);
        }
    }
    const ChVector3d origin = site.ToLocal(site_lon, site_lat, 100.0);
    EXPECT_NEAR(origin.Length(), 0.0, 1e-9);
}

TEST(ChPlanetSurface, GridMatchesPointSamples) {
    auto surface = MakeSurface();
    const int n = 5;
    const double step = 1e-4;
    std::vector<double> grid;
    surface->GetElevationGrid(site_lon, site_lat, step, step, n, grid);
    ASSERT_EQ(grid.size(), static_cast<size_t>(n * n));
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            EXPECT_NEAR(grid[j * n + i], surface->GetElevation(site_lon + i * step, site_lat + j * step), 1e-6);
}

TEST(PlanetTerrain, HeightAndNormal) {
    ChSystemNSC sys;
    auto surface = MakeSurface();
    ChSiteFrame site = surface->MakeSiteFrame(site_lon, site_lat);
    PlanetTerrain terrain(&sys, surface, site);
    terrain.Initialize();

    EXPECT_NEAR(terrain.GetHeight(ChVector3d(0, 0, 50)), 0.0, 1e-9);
    for (double x : {-7.0, 3.5}) {
        const ChVector3d loc(x, 2.0, 10.0);
        double lon, lat;
        site.ToLonLat(loc.x(), loc.y(), lon, lat);
        const double expected = surface->GetElevation(lon, lat) - site.GetOriginElevation();
        EXPECT_NEAR(terrain.GetHeight(loc), expected, 1e-9);
        EXPECT_NEAR(terrain.GetPoint(loc).z(), expected, 1e-9);
        const ChVector3d normal = terrain.GetNormal(loc);
        EXPECT_NEAR(normal.Length(), 1.0, 1e-9);
        EXPECT_GT(normal.z(), 0.5);
    }
    EXPECT_FLOAT_EQ(terrain.GetCoefficientFriction(ChVector3d(0, 0, 0)), 0.8f);
}

TEST(PlanetTerrain, PatchFollowsLocation) {
    ChSystemNSC sys;
    auto surface = MakeSurface();
    ChSiteFrame site = surface->MakeSiteFrame(site_lon, site_lat);
    PlanetTerrain terrain(&sys, surface, site);
    terrain.SetPatchSize(20.0);
    terrain.SetPatchResolution(0.5);
    terrain.SetRebuildMargin(4.0);
    terrain.Initialize();

    auto first = terrain.GetGroundBody();
    ASSERT_TRUE(first);
    EXPECT_EQ(sys.GetBodies().size(), 1u);

    EXPECT_FALSE(terrain.UpdatePatch(ChVector3d(5.9, 0, 0)));
    EXPECT_EQ(terrain.GetGroundBody(), first);

    EXPECT_TRUE(terrain.UpdatePatch(ChVector3d(6.1, 0, 0)));
    EXPECT_NE(terrain.GetGroundBody(), first);
    EXPECT_NEAR(terrain.GetPatchCenter().x(), 6.1, 1e-12);
    EXPECT_EQ(sys.GetBodies().size(), 1u);
}

TEST(PlanetTerrain, BallRestsOnPatch) {
    ChSystemNSC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, -moon::kGravity));
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    auto surface = MakeSurface();
    ChSiteFrame site = surface->MakeSiteFrame(site_lon, site_lat);
    PlanetTerrain terrain(&sys, surface, site);
    terrain.SetPatchSize(10.0);
    terrain.Initialize();

    const double radius = 0.3;
    auto material = ChContactMaterial::DefaultMaterial(sys.GetContactMethod());
    auto ball = chrono_types::make_shared<ChBodyEasySphere>(radius, 1000, true, true, material);
    ball->SetPos(ChVector3d(0, 0, 0.5));
    sys.Add(ball);
    while (sys.GetChTime() < 2.0)
        sys.DoStepDynamics(1e-3);

    const double gap = ball->GetPos().z() - terrain.GetHeight(ball->GetPos()) - radius;
    EXPECT_NEAR(gap, 0.0, 0.05);
}

TEST(PlanetSCMTerrain, FunctorMatchesSurface) {
    auto surface = MakeSurface();
    ChSiteFrame site = surface->MakeSiteFrame(site_lon, site_lat);
    PlanetSCMHeightFunctor functor(surface, site);
    const double delta = 0.05;

    EXPECT_EQ(functor.GetNumCachedNodes(), 0u);
    std::vector<ChVector2i> locs = {{0, 0}, {13, -7}, {-40, 25}};
    for (const auto& loc : locs) {
        double lon, lat;
        site.ToLonLat(loc.x() * delta, loc.y() * delta, lon, lat);
        const double expected = surface->GetElevation(lon, lat) - site.GetOriginElevation();
        EXPECT_NEAR(functor.GetInitHeight(loc, delta), expected, 1e-9);
        EXPECT_NEAR(functor.GetInitHeight(loc, delta), expected, 1e-9);  // memoized
    }
    EXPECT_EQ(functor.GetNumCachedNodes(), locs.size());

    std::vector<double> batch;
    functor.GetInitHeights(locs, delta, batch);
    ASSERT_EQ(batch.size(), locs.size());
    for (size_t i = 0; i < locs.size(); ++i)
        EXPECT_NEAR(batch[i], functor.GetInitHeight(locs[i], delta), 0.0);
}

TEST(PlanetSCMTerrain, WheelSinks) {
    ChSystemNSC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, -moon::kGravity));
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    auto surface = MakeSurface();
    ChSiteFrame site = surface->MakeSiteFrame(site_lon, site_lat);

    const double radius = 0.25;
    const double width = 0.2;
    // SCM finds the wheel by ray-casting the collision system, so the wheel needs a collision shape.
    auto wheel_mat = ChContactMaterial::DefaultMaterial(sys.GetContactMethod());
    auto wheel = chrono_types::make_shared<ChBodyEasyCylinder>(ChAxis::Y, radius, width, 2000, true, true, wheel_mat);
    wheel->SetPos(ChVector3d(0, 0, radius + 0.1));
    sys.Add(wheel);

    PlanetSCMTerrain terrain(&sys, surface, site);
    EXPECT_THROW(terrain.Initialize(PlanetSCMTerrain::Params{}, {}), std::invalid_argument);

    PlanetSCMTerrain::Params params;
    terrain.Initialize(params, {{wheel, radius, width}});
    while (sys.GetChTime() < 1.0)
        sys.DoStepDynamics(1e-3);

    EXPECT_GT(terrain.GetNumCachedNodes(), 0u);
    EXPECT_GT(terrain.GetContactForce(wheel), 0.0);
    EXPECT_LT(wheel->GetPos().z(), radius);  // below the undeformed surface: it sank

    terrain.RequestRutStats();
    std::optional<PlanetSCMTerrain::RutStats> stats;
    for (int i = 0; i < 5000 && !(stats = terrain.TakeRutStats()); ++i)
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    ASSERT_TRUE(stats.has_value());
    EXPECT_GT(stats->deformed, 0u);
    EXPECT_GT(stats->max_depth_m, 0.0);
}
