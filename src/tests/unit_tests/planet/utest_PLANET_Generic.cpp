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
// Tests of the body-agnostic Planet API: arbitrary bodies, raster options,
// fallbacks, relief layers and the quadtree, with no Moon assumptions.
//
// =============================================================================

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "chrono/core/ChTypes.h"

#include "chrono_planet/ChPlanetBody.h"
#include "chrono_planet/ChPlanetSurface.h"
#include "chrono_planet/filters/ChDeformationFilter.h"
#include "chrono_planet/lod/ChPlanetQuadtree.h"
#include "chrono_planet/lod/SphericalCoordinates.h"
#include "chrono_planet/lod/CartesianCoordinates.h"
#include "chrono_planet/lod/QuadtreeTile.h"
#include "chrono_planet/lod/TileMeshBuilder.h"
#include "chrono_planet/planets/mars/ChMars.h"
#include "chrono_planet/planets/moon/ChMoonDem.h"
#include "chrono_planet/procedural/ChBaseReliefLayer.h"
#include "chrono_planet/procedural/ChCraterLayer.h"
#include "chrono_planet/procedural/ChRockLayer.h"
#include "chrono_planet/procedural/ChRoughnessLayer.h"

#include "planet_test_utils.h"

using namespace chrono;
using namespace chrono::planet;
using planet_test::kPi;

namespace {

// A small airless body, far from the Moon's radius, so any leftover lunar constant shows.
const ChPlanetBody kBody("Testworld", 250000.0, 0.28);

double Terrain(double lon, double lat) {
    return 100.0 * std::sin(lon * kPi / 180) * std::cos(lat * kPi / 180) + 20.0 * std::sin(lat * kPi / 90);
}

// The cubic B-spline smooths the samples, so a 1-degree raster of Terrain() reads back within this.
constexpr double kReadTol = 0.1;

std::string Path(const std::string& name) {
    return planet_test::TempDir() + "/" + name;
}

// A global 1-degree raster of Terrain() on kBody, written once per name.
std::string GlobalRaster(const std::string& name = "global.tif") {
    const std::string path = Path(name);
    planet_test::WriteGlobal(path, kBody.GetRadius(), 1.0, -180.0, Terrain);
    return path;
}

}  // namespace

// -----------------------------------------------------------------------------

TEST(ChPlanetBody, Geometry) {
    EXPECT_THROW(ChPlanetBody("bad", 0.0), std::invalid_argument);
    EXPECT_NEAR(kBody.GetMetersPerDegree(), 250000.0 * kPi / 180, 1e-9);
    EXPECT_DOUBLE_EQ(kBody.GetSurfaceGravity(), 0.28);

    const ChVector3d p = kBody.ToCartesian(40.0, -25.0, 1200.0);
    EXPECT_NEAR(p.Length(), 251200.0, 1e-6);
    double lon, lat, h;
    kBody.ToGeographic(p, lon, lat, h);
    EXPECT_NEAR(lon, 40.0, 1e-9);
    EXPECT_NEAR(lat, -25.0, 1e-9);
    EXPECT_NEAR(h, 1200.0, 1e-6);
}

TEST(ChSiteFrame, UsesBodyRadius) {
    const ChSiteFrame site(kBody, 10.0, 5.0, 0.0);
    EXPECT_DOUBLE_EQ(site.GetRadius(), kBody.GetRadius());
    // One degree north is one degree of this body's meridian.
    EXPECT_NEAR(site.ToLocal(10.0, 6.0, 0.0).y(), kBody.GetMetersPerDegree(), 1e-6);
}

// -----------------------------------------------------------------------------

TEST(ChPlanetSurface, FlatByDefault) {
    ChPlanetSurface surface(kBody);
    EXPECT_EQ(surface.GetElevation(12.0, 34.0), 0.0);
    std::vector<double> h;
    surface.GetElevationGrid(0.0, 0.0, 0.01, 0.01, 5, h);
    for (double v : h)
        EXPECT_EQ(v, 0.0);
}

TEST(ChPlanetSurface, ReadsAnyBodyDem) {
    ChPlanetSurface surface(kBody);
    surface.AddGeoTiff({GlobalRaster(), 0, 30});
    for (double lon : {-150.0, -30.5, 0.0, 77.25, 170.0})
        for (double lat : {-60.0, -5.5, 0.0, 33.3, 70.0})
            EXPECT_NEAR(surface.GetElevation(lon, lat), Terrain(lon, lat), kReadTol) << lon << " " << lat;
    EXPECT_EQ(surface.GetGeoTiffStack()->GetNumSources(), 1);
}

TEST(ChPlanetSurface, CustomElevationSource) {
    class Ramp : public ChElevationSampler {
      public:
        std::optional<double> GetHeight(double lon, double lat, int) const override {
            if (lon < 0)
                return std::nullopt;
            return 10.0 * lon + lat;
        }
    };
    ChPlanetSurface surface(kBody);
    surface.SetSampler(chrono_types::make_shared<Ramp>());
    surface.SetFallback(chrono_types::make_shared<ChElevationFunction>([](double, double) { return -5.0; }));
    EXPECT_DOUBLE_EQ(surface.GetElevation(3.0, 2.0), 32.0);
    EXPECT_DOUBLE_EQ(surface.GetElevation(-3.0, 2.0), -5.0);
    std::vector<double> h;
    surface.GetElevationGrid(-0.5, 1.0, 0.25, 0.5, 4, h);  // columns at -0.5, -0.25, 0, 0.25
    EXPECT_DOUBLE_EQ(h[0], -5.0);
    EXPECT_DOUBLE_EQ(h[1], -5.0);
    EXPECT_DOUBLE_EQ(h[2], 1.0);
    EXPECT_DOUBLE_EQ(h[3], 3.5);
    // AddGeoTiff cannot mix with a different kind of source.
    EXPECT_THROW(surface.AddGeoTiff({GlobalRaster(), 0, 30}), std::logic_error);
}

TEST(ChPlanetSurface, MissingFileThrows) {
    ChPlanetSurface surface(kBody);
    try {
        surface.AddGeoTiff({Path("does_not_exist.tif"), 0, 30});
        FAIL() << "no exception";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("does_not_exist.tif"), std::string::npos);
    }
}

TEST(ChPlanetSurface, SiteFrameOnSurface) {
    ChPlanetSurface surface(kBody);
    surface.AddGeoTiff({GlobalRaster(), 0, 30});
    const ChSiteFrame site = surface.MakeSiteFrame(20.0, 10.0);
    EXPECT_DOUBLE_EQ(site.GetOriginElevation(), surface.GetElevation(20.0, 10.0));
    EXPECT_DOUBLE_EQ(site.GetRadius(), kBody.GetRadius());
}

// -----------------------------------------------------------------------------

TEST(ChGeoTiffSource, ValueConversions) {
    const auto write = [](const std::string& name, auto stored, planet_test::RasterOptions opt = {}) {
        const std::string path = Path(name);
        planet_test::WriteGlobal(path, kBody.GetRadius(), 1.0, -180.0,
                                 [&](double lon, double lat) { return stored(Terrain(lon, lat)); }, opt);
        return path;
    };
    const double lon = 42.5, lat = -12.25, expected = Terrain(lon, lat);

    {  // kilometers, and an offset datum
        ChPlanetSurface s(kBody);
        ChGeoTiffSource src{write("km.tif", [](double h) { return (h - 50.0) / 1000.0; }), 0, 30};
        src.scale = 1000.0;
        src.offset = 50.0;
        s.AddGeoTiff(src);
        EXPECT_NEAR(s.GetElevation(lon, lat), expected, kReadTol);
    }
    {  // distances from the body center
        ChPlanetSurface s(kBody);
        ChGeoTiffSource src{write("radii.tif", [](double h) { return kBody.GetRadius() + h; }), 0, 30};
        src.values_are_radii = true;
        s.AddGeoTiff(src);
        // Float32 radii carry ~2 cm of precision at this radius.
        EXPECT_NEAR(s.GetElevation(lon, lat), expected, kReadTol + 0.05);
    }
    {  // the band's own scale/offset metadata, applied by default and skippable
        planet_test::RasterOptions opt;
        opt.band_scale = 0.5;
        opt.band_offset = 10.0;
        const std::string path = write("bandscale.tif", [](double h) { return (h - 10.0) / 0.5; }, opt);
        ChPlanetSurface s(kBody);
        s.AddGeoTiff({path, 0, 30});
        EXPECT_NEAR(s.GetElevation(lon, lat), expected, kReadTol);

        ChPlanetSurface raw(kBody);
        ChGeoTiffSource src{path, 0, 30};
        src.use_band_scale_offset = false;
        raw.AddGeoTiff(src);
        EXPECT_NEAR(raw.GetElevation(lon, lat), (expected - 10.0) / 0.5, 2 * kReadTol);
    }
    {  // data in a later band
        planet_test::RasterOptions opt;
        opt.bands = 2;
        const std::string path = write("band2.tif", [](double h) { return h; }, opt);
        ChPlanetSurface s(kBody);
        ChGeoTiffSource src{path, 0, 30};
        src.band = 2;
        s.AddGeoTiff(src);
        EXPECT_NEAR(s.GetElevation(lon, lat), expected, kReadTol);
        ChGeoTiffSource bad{path, 0, 30};
        bad.band = 3;
        EXPECT_THROW(s.AddGeoTiff(bad), std::runtime_error);
    }
}

TEST(ChGeoTiffSource, LongitudeConventionsAndMissingSrs) {
    // 0..360 longitudes, and a raster with no coordinate system at all.
    const std::string east = Path("east360.tif"), bare = Path("nosrs.tif");
    planet_test::WriteGlobal(east, kBody.GetRadius(), 1.0, 0.0, Terrain);
    planet_test::WriteGlobal(bare, kBody.GetRadius(), 1.0, -180.0, Terrain, {}, false);
    for (const std::string& path : {east, bare}) {
        ChPlanetSurface s(kBody);
        s.AddGeoTiff({path, 0, 30});
        // Avoid the last pixel before each raster's wrap edge: global rasters are not bridged across it.
        for (double lon : {-100.0, -1.5, 45.0, 135.0})
            EXPECT_NEAR(s.GetElevation(lon, 15.0), Terrain(lon, 15.0), kReadTol) << path << " " << lon;
    }
}

TEST(ChGeoTiffSource, ProjectedRasterOnAnyBody) {
    // An equirectangular raster in meters on kBody's sphere, around (60, 30).
    OGRSpatialReference proj;
    proj.SetProjCS("Test_Equirectangular");
    proj.SetGeogCS("GCS_Test", "D_Test", "Test_Sphere", kBody.GetRadius(), 0.0);
    proj.SetEquirectangular2(0.0, 0.0, 0.0, 0.0, 0.0);
    const double R = kBody.GetRadius(), px = 200.0;
    const int n = 200;
    const double cx = R * 60.0 * kPi / 180, cy = R * 30.0 * kPi / 180;
    const double gt[6] = {cx - 0.5 * n * px, px, 0.0, cy + 0.5 * n * px, 0.0, -px};
    std::vector<float> data(n * n);
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c) {
            // The stack places a pixel's value at its geotransform corner (see ChGeoTiffSource).
            const double x = gt[0] + c * px, y = gt[3] - r * px;
            data[r * n + c] = static_cast<float>(Terrain(x / R * 180 / kPi, y / R * 180 / kPi));
        }
    const std::string path = Path("projected.tif");
    planet_test::WriteRaster(path, n, n, gt, &proj, data);

    ChPlanetSurface s(kBody);
    s.AddGeoTiff({path, 0, 30});
    EXPECT_NEAR(s.GetElevation(60.2, 29.9), Terrain(60.2, 29.9), 0.01);
    EXPECT_FALSE(s.GetDataElevation(70.0, 30.0, 15).has_value());  // outside the raster
}

TEST(ChGeoTiffSource, BaseHolesFallBack) {
    // A global base with a nodata cap north of 60 degrees: the cap reads as the fallback, not the sentinel.
    planet_test::RasterOptions opt;
    opt.nodata = -32768.0;
    const std::string path = Path("holes.tif");
    planet_test::WriteGlobal(path, kBody.GetRadius(), 1.0, -180.0,
                             [](double lon, double lat) { return lat > 60.0 ? -32768.0 : Terrain(lon, lat); }, opt);
    ChPlanetSurface s(kBody);
    s.AddGeoTiff({path, 0, 30});
    s.SetFallback(chrono_types::make_shared<ChElevationFunction>([](double, double) { return -777.0; }));
    EXPECT_DOUBLE_EQ(s.GetElevation(10.0, 75.0), -777.0);
    EXPECT_NEAR(s.GetElevation(10.0, 20.0), Terrain(10.0, 20.0), kReadTol);
    std::vector<double> h;
    s.GetElevationGrid(10.0, 74.0, 0.5, 0.5, 4, h);
    for (double v : h)
        EXPECT_DOUBLE_EQ(v, -777.0);

    // Overriding nodata makes the sentinel an ordinary value.
    ChPlanetSurface raw(kBody);
    ChGeoTiffSource src{path, 0, 30};
    src.nodata = 12345.0;
    src.min_valid = -1e9;
    raw.AddGeoTiff(src);
    EXPECT_NEAR(raw.GetElevation(10.0, 75.0), -32768.0, 1.0);
}

TEST(ChGeoTiffSource, OverlayByZoom) {
    // A coarse global base and a regional overlay offset by +100 m, used from zoom 5 on.
    const std::string base = GlobalRaster("zoom_base.tif");
    const std::string overlay = Path("zoom_overlay.tif");
    {
        const double gt[6] = {10.0, 0.01, 0.0, 21.0, 0.0, -0.01};
        std::vector<float> data(200 * 200);
        for (int r = 0; r < 200; ++r)
            for (int c = 0; c < 200; ++c)
                data[r * 200 + c] = static_cast<float>(Terrain(10.0 + c * 0.01, 21.0 - r * 0.01) + 100.0);
        const OGRSpatialReference geo = planet_test::SphereGeog(kBody.GetRadius());
        planet_test::WriteRaster(overlay, 200, 200, gt, &geo, data);
    }
    ChPlanetSurface s(kBody);
    s.AddGeoTiff({base, 0, 4});
    s.AddGeoTiff({overlay, 5, 30});
    EXPECT_NEAR(s.GetElevationAtZoom(11.0, 20.0, 3), Terrain(11.0, 20.0), kReadTol);
    EXPECT_NEAR(s.GetElevationAtZoom(11.0, 20.0, 10), Terrain(11.0, 20.0) + 100.0, kReadTol);
    // Outside the overlay the base answers at every zoom.
    EXPECT_NEAR(s.GetElevationAtZoom(-50.0, 20.0, 10), Terrain(-50.0, 20.0), kReadTol);
}

// -----------------------------------------------------------------------------

TEST(ChReliefLayer, DefaultsAreFlat) {
    // No body-specific numbers hide in the defaults: every layer is zero until configured.
    const ChBaseReliefLayer base(ChBaseReliefLayer::Params{});
    const ChCraterLayer craters(kBody, ChCraterLayer::Params{});
    const ChRockLayer rocks(kBody, ChRockLayer::Params{});
    const ChRoughnessLayer roughness(kBody, ChRoughnessLayer::Params{});
    for (double lon : {-20.0, 0.001, 33.0})
        for (double spacing : {1e-3, 1e-7}) {
            EXPECT_EQ(base.GetHeight(lon, 10.0, spacing), 0.0);
            EXPECT_EQ(craters.GetHeight(lon, 10.0, spacing), 0.0);
            EXPECT_EQ(rocks.GetHeight(lon, 10.0, spacing), 0.0);
            EXPECT_EQ(roughness.GetHeight(lon, 10.0, spacing), 0.0);
        }
    EXPECT_TRUE(rocks.Query(0, 0, 0.01, 0.01, 0).empty());
    EXPECT_EQ(rocks.GetMaxRadius(), 0.0);
}

TEST(ChReliefLayer, InvalidParams) {
    ChCraterLayer::Params c;
    c.diameters_km = {1.0, 2.0};
    c.densities = {0.1, 0.1};
    EXPECT_THROW(ChCraterLayer(kBody, c), std::invalid_argument);  // not decreasing
    c.diameters_km = {2.0, 1.0};
    c.densities = {0.1};
    EXPECT_THROW(ChCraterLayer(kBody, c), std::invalid_argument);  // length mismatch
    ChRockLayer::Params r;
    r.coverage = 1.5;
    EXPECT_THROW(ChRockLayer(kBody, r), std::invalid_argument);
    ChRoughnessLayer::Params g;
    g.clod_variance_share = 2;
    EXPECT_THROW(ChRoughnessLayer(kBody, g), std::invalid_argument);
}

TEST(ChReliefLayer, GridMatchesPointsOnAnyBody) {
    // Lunar-style parameters on a different body: grid and point evaluation must agree.
    auto surface = chrono_types::make_shared<ChPlanetSurface>(kBody, 12);
    surface->AddLayer(chrono_types::make_shared<ChBaseReliefLayer>(moon::BaseReliefParams()));
    surface->AddLayer(chrono_types::make_shared<ChCraterLayer>(kBody, moon::CraterParams()));
    surface->AddLayer(chrono_types::make_shared<ChRockLayer>(kBody, moon::RockParams()));
    surface->AddLayer(chrono_types::make_shared<ChRoughnessLayer>(kBody, moon::RoughnessParams()));

    const double spacing = surface->GetSampleSpacing();
    for (double lon0 : {45.0, 179.9995}) {  // the second grid crosses the dateline
        const ChGeoGrid grid{lon0, -8.0, spacing, spacing, 17, spacing};
        std::vector<double> h;
        surface->GetElevationGrid(grid, 12, h);
        double max_relief = 0;
        for (int j = 0; j < grid.n; ++j)
            for (int i = 0; i < grid.n; ++i) {
                const double lon = lon0 + i * grid.step_lon, lat = grid.lat0 + j * grid.step_lat;
                const double point = surface->GetElevationAtZoom(lon, lat, 12, spacing);
                EXPECT_NEAR(h[j * grid.n + i], point, 1e-6) << lon0 << " " << i << " " << j;
                max_relief = std::max(max_relief, std::abs(point));
            }
        EXPECT_GT(max_relief, 0.0);
    }
}

TEST(ChReliefLayer, SeedChangesRealization) {
    ChCraterLayer::Params p = moon::CraterParams();
    const ChCraterLayer a(kBody, p);
    p.seed = 7;
    const ChCraterLayer b(kBody, p);
    double diff = 0;
    for (int i = 0; i < 50; ++i)
        diff += std::abs(a.GetHeight(0.1 * i, 3.0, 1e-4) - b.GetHeight(0.1 * i, 3.0, 1e-4));
    EXPECT_GT(diff, 0.0);
}

TEST(ChReliefLayer, CraterSizesAreInMeters) {
    // The same crater parameters on a body twice the size give craters of the same physical size,
    // so the relief as a function of ground distance is unchanged in distribution: at equal ground
    // spacing, both resolve the same number of classes.
    const ChPlanetBody big("Big", 2 * kBody.GetRadius());
    const ChCraterLayer small_layer(kBody, moon::CraterParams());
    const ChCraterLayer big_layer(big, moon::CraterParams());
    // 10 m ground spacing on each body.
    const double small_spacing = 10.0 / kBody.GetMetersPerDegree(), big_spacing = 10.0 / big.GetMetersPerDegree();
    double small_rms = 0, big_rms = 0;
    const int n = 400;
    for (int i = 0; i < n; ++i) {
        const double hs = small_layer.GetHeight(i * 50 * small_spacing, 0.5, small_spacing);
        const double hb = big_layer.GetHeight(i * 50 * big_spacing, 0.5, big_spacing);
        small_rms += hs * hs;
        big_rms += hb * hb;
    }
    // Different realizations, similar statistics.
    EXPECT_GT(small_rms, 0.0);
    EXPECT_GT(big_rms, 0.0);
    EXPECT_LT(std::abs(std::log(small_rms / big_rms)), std::log(10.0));
}

TEST(ChPlanetSurface, FindFilter) {
    auto surface = moon::CreateSurface(15);
    ASSERT_TRUE(surface->FindFilter<ChRockLayer>());
    EXPECT_GT(surface->FindFilter<ChRockLayer>()->GetNumClasses(), 0);
    ChPlanetSurface flat(kBody);
    EXPECT_FALSE(flat.FindFilter<ChCraterLayer>());
}

// -----------------------------------------------------------------------------

TEST(Presets, BodiesAreDistinct) {
    EXPECT_DOUBLE_EQ(moon::Body().GetRadius(), moon::kRadius);
    EXPECT_DOUBLE_EQ(mars::Body().GetRadius(), mars::kRadius);
    auto surface = mars::CreateSurface(std::vector<ChGeoTiffSource>(), 12);
    EXPECT_EQ(surface->GetBody().GetName(), "Mars");
    EXPECT_NE(surface->GetElevation(137.4, -4.6), 0.0);  // relief from the crater and rock layers
}

TEST(ChPlanetQuadtree, AnyBody) {
    auto surface = chrono_types::make_shared<ChPlanetSurface>(kBody, 12);
    surface->AddGeoTiff({GlobalRaster("quadtree.tif"), 0, 30});
    surface->SetRootTileSize(8.0);
    ChPlanetQuadtree quadtree(surface, 1);
    const ChVector3d cam = kBody.ToCartesian(25.0, 15.0, 2000.0);
    quadtree.Update(cam, 0.0);
    EXPECT_GT(quadtree.GetMeshes().size(), 0u);
    EXPECT_NEAR(quadtree.GetCameraElevation(), 2000.0, 1e-6);
    for (int zoom : {0, 6, 12}) {
        const auto h = quadtree.GetElevation(25.1, 14.9, zoom);
        ASSERT_TRUE(h.has_value());
        EXPECT_DOUBLE_EQ(*h, surface->GetElevationAtZoom(25.1, 14.9, zoom));
    }
    // Meshes sit on this body's sphere: every vertex within the terrain's range of the radius.
    for (const ChTileMesh* mesh : quadtree.GetMeshes()) {
        const ChTileMesh& m = *mesh;
        EXPECT_GT(m.minElevation, -200.0);
        EXPECT_LT(m.maxElevation, 200.0);
    }
}

// -----------------------------------------------------------------------------
// Filter chains

TEST(ChFilterChain, OrderInsertRemoveFind) {
    auto chain = chrono_types::make_shared<ChFilterChain>();
    auto scale = chrono_types::make_shared<ChScaleFilter>(2.0);
    auto shift = chrono_types::make_shared<ChScaleFilter>(1.0, 10.0);
    chain->AddFilter(scale);
    chain->AddFilter(shift);
    EXPECT_DOUBLE_EQ(chain->Apply(0, 0, 1e-3, 5.0), 20.0);  // (5 * 2) + 10
    chain->RemoveFilter(shift);
    chain->InsertFilter(0, shift);
    EXPECT_DOUBLE_EQ(chain->Apply(0, 0, 1e-3, 5.0), 30.0);  // (5 + 10) * 2
    EXPECT_FALSE(chain->RemoveFilter(chrono_types::make_shared<ChScaleFilter>(1.0)));
    EXPECT_THROW(chain->AddFilter(chain), std::invalid_argument);

    // Nested chains are searched.
    auto outer = chrono_types::make_shared<ChFilterChain>();
    outer->AddFilter(chain);
    outer->AddFilter(chrono_types::make_shared<ChClampFilter>(-1.0, 1.0));
    EXPECT_EQ(outer->Find<ChScaleFilter>(), shift);
    EXPECT_TRUE(outer->Find<ChClampFilter>());
    EXPECT_FALSE(outer->Find<ChRegionFilter>());
}

TEST(ChFilterChain, PresetChainIsExtensible) {
    // The Moon's chain with a user filter appended: vertical exaggeration of the whole terrain.
    auto plain = moon::CreateSurface(15);
    auto tall = moon::CreateSurface(15);
    tall->AddFilter(chrono_types::make_shared<ChScaleFilter>(3.0));
    for (double lon : {10.0, 30.75})
        EXPECT_NEAR(tall->GetElevation(lon, 20.0), 3.0 * plain->GetElevation(lon, 20.0), 1e-9);
    EXPECT_EQ(moon::FilterChain()->GetNumFilters(), 4u);
    EXPECT_EQ(mars::FilterChain()->GetNumFilters(), 2u);
    EXPECT_NE(moon::FilterChain(), moon::FilterChain());  // independent copies
}

TEST(ChFilterChain, UserFilterAndGridMatchesPoints) {
    // A user filter that terraces the ground in 25 m steps, after the relief.
    class Terraces : public ChSurfaceFilter {
      public:
        double Apply(double, double, double, double h) const override { return 25.0 * std::floor(h / 25.0); }
    };
    ChPlanetSurface surface(kBody, 10);
    surface.AddGeoTiff({GlobalRaster("filters.tif"), 0, 30});
    surface.AddLayer(chrono_types::make_shared<ChCraterLayer>(kBody, moon::CraterParams()));
    surface.AddFilter(chrono_types::make_shared<ChFunctionFilter>(
        [](double lon, double, double, double h) { return h + 0.5 * lon; }));
    surface.AddFilter(chrono_types::make_shared<Terraces>());
    surface.AddFilter(chrono_types::make_shared<ChClampFilter>(-60.0, 60.0));

    for (double lon0 : {20.0, 179.9}) {  // the second grid crosses the dateline
        const ChGeoGrid grid{lon0, 5.0, 0.05, 0.05, 9, surface.GetSampleSpacing()};
        std::vector<double> h;
        surface.GetElevationGrid(grid, 10, h);
        for (int j = 0; j < grid.n; ++j)
            for (int i = 0; i < grid.n; ++i) {
                const double v = h[j * grid.n + i];
                EXPECT_EQ(std::fmod(v, 25.0) == 0.0 || v == 60.0 || v == -60.0, true) << v;
                EXPECT_DOUBLE_EQ(v, surface.GetElevationAtZoom(lon0 + i * 0.05, 5.0 + j * 0.05, 10, grid.spacing));
            }
    }
}

TEST(ChRegionFilter, WeightsAndBlend) {
    auto flatten = chrono_types::make_shared<ChFunctionFilter>([](double, double, double, double) { return 0.0; });
    // A landing pad across the dateline: longitudes 179..-179, latitudes 0..1, feathered over 0.5 degrees.
    const ChRegionFilter pad(flatten, 179.0, 0.0, -179.0, 1.0, 0.5);
    EXPECT_DOUBLE_EQ(pad.GetWeight(179.5, 0.5), 1.0);
    EXPECT_DOUBLE_EQ(pad.GetWeight(-179.5, 0.5), 1.0);
    EXPECT_DOUBLE_EQ(pad.GetWeight(178.0, 0.5), 0.0);
    const double w = pad.GetWeight(178.75, 0.5);
    EXPECT_GT(w, 0.0);
    EXPECT_LT(w, 1.0);
    EXPECT_DOUBLE_EQ(pad.Apply(179.5, 0.5, 1e-3, 100.0), 0.0);
    EXPECT_DOUBLE_EQ(pad.Apply(170.0, 0.5, 1e-3, 100.0), 100.0);
    EXPECT_NEAR(pad.Apply(178.75, 0.5, 1e-3, 100.0), 100.0 * (1 - w), 1e-12);

    ChPlanetSurface surface(kBody);
    surface.AddGeoTiff({GlobalRaster("region.tif"), 0, 30});
    surface.AddFilter(chrono_types::make_shared<ChRegionFilter>(flatten, 10.0, 10.0, 11.0, 11.0, 0.25));
    const ChGeoGrid grid{9.5, 9.5, 0.1, 0.1, 21, 1e-3};
    std::vector<double> h;
    surface.GetElevationGrid(grid, 15, h);
    for (int j = 0; j < grid.n; ++j)
        for (int i = 0; i < grid.n; ++i)
            EXPECT_NEAR(h[j * grid.n + i], surface.GetElevationAtZoom(9.5 + i * 0.1, 9.5 + j * 0.1, 15, 1e-3), 1e-9);
    EXPECT_DOUBLE_EQ(surface.GetElevation(10.5, 10.5), 0.0);
}

// -----------------------------------------------------------------------------
// Sampler chains

TEST(ChSamplerChain, FirstWithDataWins) {
    auto east = chrono_types::make_shared<ChElevationFunction>([](double, double) { return 1.0; });
    class WestOnly : public ChElevationSampler {
      public:
        std::optional<double> GetHeight(double lon, double, int) const override {
            return lon < 0 ? std::optional<double>(-1.0) : std::nullopt;
        }
    };
    auto chain = chrono_types::make_shared<ChSamplerChain>();
    chain->AddSampler(chrono_types::make_shared<WestOnly>());
    chain->AddSampler(east);
    EXPECT_EQ(chain->GetHeight(-10, 0, 0), -1.0);
    EXPECT_EQ(chain->GetHeight(10, 0, 0), 1.0);

    ChPlanetSurface surface(kBody);
    surface.SetSampler(chain);
    std::vector<double> h;
    surface.GetElevationGrid(-0.2, 0.0, 0.1, 0.1, 5, h);
    EXPECT_DOUBLE_EQ(h[0], -1.0);
    EXPECT_DOUBLE_EQ(h[1], -1.0);
    EXPECT_DOUBLE_EQ(h[2], 1.0);
    EXPECT_DOUBLE_EQ(h[4], 1.0);

    // A GeoTIFF stack in a chain, in front of an analytic datum.
    auto stack = chrono_types::make_shared<ChGeoTiffStack>(kBody);
    stack->AddSource({Path("projected.tif"), 0, 30});  // written by ProjectedRasterOnAnyBody, covers ~(60, 30)
    ChPlanetSurface mixed(kBody);
    mixed.SetSampler(chrono_types::make_shared<ChSamplerChain>(
        std::vector<std::shared_ptr<ChElevationSampler>>{stack, chrono_types::make_shared<ChElevationFunction>(
                                                                    [](double, double) { return -500.0; })}));
    EXPECT_NEAR(mixed.GetElevation(60.2, 29.9), Terrain(60.2, 29.9), 0.01);
    EXPECT_DOUBLE_EQ(mixed.GetElevation(-60.0, 0.0), -500.0);
}

// -----------------------------------------------------------------------------
// Moon DEM resources

TEST(MoonDem, Resources) {
    const ChGeoTiffSource apollo = moon::GetDem(moon::Dem::APOLLO17_LANDING_SITE, "/some/dir");
    EXPECT_EQ(apollo.path, "/some/dir/ldem_1024_apollo_region.tif");
    EXPECT_EQ(apollo.min_zoom, 5);
    EXPECT_EQ(moon::GetDem(moon::Dem::GLOBAL_LOW_RES, "/d/").path, "/d/ldem_4_global.tif");
    EXPECT_FALSE(moon::GetDemDescription(moon::Dem::GLOBAL).empty());

    ChPlanetSurface surface(moon::Body());
    EXPECT_FALSE(moon::IsDemAvailable(moon::Dem::GLOBAL, "/nonexistent"));
    EXPECT_THROW(moon::AddDem(surface, moon::Dem::GLOBAL, "/nonexistent"), std::runtime_error);
}

TEST(MoonDem, ShippedResourcesLoad) {
    if (!moon::IsDemAvailable(moon::Dem::GLOBAL_LOW_RES) || !moon::IsDemAvailable(moon::Dem::APOLLO17_LANDING_SITE))
        GTEST_SKIP() << "Moon DEM resources not found under " << moon::DataDir();
    auto surface = moon::CreateSurface({moon::Dem::GLOBAL_LOW_RES, moon::Dem::APOLLO17_LANDING_SITE}, 15);
    EXPECT_EQ(surface->GetGeoTiffStack()->GetNumSources(), 2);
    // At the Apollo 17 site the regional DEM answers: the value the module produced with it before bodies
    // became configurable.
    EXPECT_NEAR(surface->GetElevation(30.75, 20.19), -2225.4056142171785, 1e-6);
    // Far from the site, the global model answers at every zoom.
    const auto far = surface->GetDataElevation(-100.0, -40.0, 12);
    ASSERT_TRUE(far.has_value());
    EXPECT_GT(*far, -10000.0);
    EXPECT_LT(*far, 10000.0);
}

// -----------------------------------------------------------------------------
// Run-time deformation

TEST(ChDeformationFilter, BilinearFadeAndChanges) {
    const ChSiteFrame site(kBody, 20.0, 10.0, 0.0);
    ChDeformationFilter filter(site, 0.1);
    const double fine = 0.05 / kBody.GetMetersPerDegree();   // 5 cm relief spacing: fully visible
    const double coarse = 1.0 / kBody.GetMetersPerDegree();  // 1 m: 10 node spacings, invisible

    EXPECT_EQ(filter.GetChanges(0, *std::make_unique<ChGeoRegion>()), 0u);
    EXPECT_EQ(filter.SetDeltas({{ChVector2i(0, 0), -0.2}, {ChVector2i(1, 0), -0.1}}), 2u);
    EXPECT_EQ(filter.SetDeltas({{ChVector2i(0, 0), -0.2005}}), 0u);  // within the 1 mm tolerance
    EXPECT_EQ(filter.GetNumNodes(), 2u);

    double lon, lat;
    site.ToLonLat(0.05, 0.0, lon, lat);  // halfway between the two nodes
    EXPECT_NEAR(filter.Apply(lon, lat, fine, 100.0), 100.0 - 0.15, 1e-9);
    EXPECT_DOUBLE_EQ(filter.Apply(lon, lat, coarse, 100.0), 100.0);
    site.ToLonLat(5.0, 5.0, lon, lat);  // far from any node
    EXPECT_DOUBLE_EQ(filter.Apply(lon, lat, fine, 100.0), 100.0);

    // Changes are reported once per version, over the touched area.
    ChGeoRegion region;
    const std::uint64_t v1 = filter.GetChanges(0, region);
    EXPECT_GT(v1, 0u);
    ASSERT_FALSE(region.IsEmpty());
    site.ToLonLat(0.05, 0.0, lon, lat);
    EXPECT_LE(region.min_lon, lon);
    EXPECT_GE(region.max_lon, lon);
    EXPECT_NEAR(region.max_spacing, 0.4 / kBody.GetMetersPerDegree(), 1e-12);
    ChGeoRegion none;
    EXPECT_EQ(filter.GetChanges(v1, none), v1);
    EXPECT_TRUE(none.IsEmpty());
    filter.SetDeltas({{ChVector2i(50, 50), -0.3}});
    ChGeoRegion later;
    EXPECT_GT(filter.GetChanges(v1, later), v1);
    site.ToLonLat(5.0, 5.0, lon, lat);
    EXPECT_LE(later.min_lon, lon);
    EXPECT_GE(later.max_lon, lon);

    // Grid and point evaluation agree.
    std::vector<double> h(9 * 9, 1.0);
    site.ToLonLat(-0.2, -0.2, lon, lat);
    const double step = 0.05 / kBody.GetMetersPerDegree();
    const ChGeoGrid grid{lon, lat, step, step, 9, fine};
    filter.ApplyGrid(grid, h);
    for (int j = 0; j < 9; ++j)
        for (int i = 0; i < 9; ++i)
            EXPECT_NEAR(h[j * 9 + i], filter.Apply(lon + i * step, lat + j * step, fine, 1.0), 1e-12);
}

TEST(ChPlanetSurface, ViewAddsFiltersWithoutTouchingPhysics) {
    auto physics = chrono_types::make_shared<ChPlanetSurface>(kBody, 17);
    physics->AddGeoTiff({GlobalRaster("view.tif"), 0, 30});
    const ChSiteFrame site = physics->MakeSiteFrame(20.0, 10.0);
    auto ruts = chrono_types::make_shared<ChDeformationFilter>(site, 0.05);
    auto drawn = physics->CreateView(ruts);

    double lon, lat;
    site.ToLonLat(1.0, 2.0, lon, lat);
    const double before = physics->GetElevation(lon, lat);
    EXPECT_DOUBLE_EQ(drawn->GetElevation(lon, lat), before);
    ruts->SetDeltas({{ChVector2i(20, 40), -0.08}});
    EXPECT_DOUBLE_EQ(physics->GetElevation(lon, lat), before);  // physics never sees the ruts
    EXPECT_NEAR(drawn->GetElevation(lon, lat), before - 0.08, 1e-9);

    // Edits to the physics chain show through the view.
    physics->AddFilter(chrono_types::make_shared<ChScaleFilter>(1.0, 5.0));
    EXPECT_NEAR(drawn->GetElevation(lon, lat), before + 5.0 - 0.08, 1e-9);
}

TEST(ChPlanetQuadtree, RebuildsTilesWhenFiltersChange) {
    auto physics = chrono_types::make_shared<ChPlanetSurface>(kBody, 17);
    physics->AddGeoTiff({GlobalRaster("rebuild.tif"), 0, 30});
    const ChSiteFrame site = physics->MakeSiteFrame(20.0, 10.0);
    auto ruts = chrono_types::make_shared<ChDeformationFilter>(site, 0.05);
    auto drawn = physics->CreateView(ruts);
    ChPlanetQuadtree quadtree(drawn, 1);

    // Settle the tree with the camera 2 m above the site.
    const ChVector3d cam = kBody.ToCartesian(20.0, 10.0, site.GetOriginElevation() + 2.0);
    double t = 0;
    for (int k = 0; k < 30; ++k)
        quadtree.Update(cam, t += ChPlanetQuadtree::kTickS);
    EXPECT_EQ(quadtree.GetNumRebuiltTiles(), 0);
    std::vector<std::uint64_t> ids_before;
    for (const ChTileMesh* m : quadtree.GetMeshes())
        ids_before.push_back(m->id);

    // A rut at the site: the next update rebuilds only fine tiles, and they carry it.
    std::vector<std::pair<ChVector2i, double>> rut;
    for (int i = -4; i <= 4; ++i)
        rut.push_back({ChVector2i(i, 0), -0.1});
    ruts->SetDeltas(rut);
    const unsigned long long version = quadtree.GetMeshSetVersion();
    quadtree.Update(cam, t += ChPlanetQuadtree::kTickS);
    EXPECT_GT(quadtree.GetNumRebuiltTiles(), 0);
    EXPECT_NE(quadtree.GetMeshSetVersion(), version);
    int replaced = 0;
    for (const ChTileMesh* m : quadtree.GetMeshes()) {
        if (std::find(ids_before.begin(), ids_before.end(), m->id) == ids_before.end()) {
            ++replaced;
            EXPECT_LT(drawn->GetSampleSpacingAtZoom(m->level), 0.2 / kBody.GetMetersPerDegree());
        }
    }
    EXPECT_GT(replaced, 0);
    double lon, lat;
    site.ToLonLat(0.0, 0.0, lon, lat);
    EXPECT_NEAR(*quadtree.GetElevation(lon, lat, 17), physics->GetElevationAtZoom(lon, lat, 17) - 0.1, 1e-9);
}

TEST(ChPlanetSurface, StaticAndDynamicStagesMatchFullEvaluation) {
    // Relief, then the deformation, then a filter after it: the split must keep every filter in order.
    auto physics = moon::CreateSurface(17);
    const ChSiteFrame site = physics->MakeSiteFrame(30.75, 20.19);
    auto ruts = chrono_types::make_shared<ChDeformationFilter>(site, 0.05);
    auto chain = chrono_types::make_shared<ChFilterChain>();
    chain->AddFilter(ruts);
    chain->AddFilter(chrono_types::make_shared<ChScaleFilter>(1.5, 2.0));
    auto drawn = physics->CreateView(chain);
    EXPECT_TRUE(drawn->HasDynamicFilters());
    EXPECT_FALSE(physics->HasDynamicFilters());
    std::vector<std::pair<ChVector2i, double>> rut;
    for (int i = -10; i <= 10; ++i)
        rut.push_back({ChVector2i(i, i / 3), -0.05});
    ruts->SetDeltas(rut);

    double lon, lat;
    site.ToLonLat(-0.5, -0.5, lon, lat);
    const double step = 0.02 / moon::Body().GetMetersPerDegree();
    const ChGeoGrid grid{lon, lat, step, step, 33, step};
    std::vector<double> full, staged;
    drawn->GetElevationGrid(grid, 17, full);
    drawn->GetStaticElevationGrid(grid, 17, staged);
    drawn->ApplyDynamicFilters(grid, staged);
    EXPECT_EQ(staged, full);
    for (int j = 0; j < grid.n; ++j)
        for (int i = 0; i < grid.n; ++i)
            EXPECT_NEAR(full[static_cast<size_t>(j) * grid.n + i],
                        drawn->GetElevationAtZoom(lon + i * step, lat + j * step, 17, step), 1e-6);
}

TEST(TileMeshBuilder, CachedRebuildMatchesFullBuild) {
    auto physics = moon::CreateSurface(17);
    const ChSiteFrame site = physics->MakeSiteFrame(30.75, 20.19);
    auto ruts = chrono_types::make_shared<ChDeformationFilter>(site, 0.05);
    auto drawn = physics->CreateView(ruts);

    // A level-17 tile over the site.
    const double size = 16.0 / (1 << 17);
    const Spherical::Boundary b{(std::floor(30.75 / size) + 0.5) * size, (std::floor((20.19 + 90) / size) + 0.5) * size - 90,
                                0.5 * size, 0.5 * size};
    const TileMeshBuilder<Spherical> builder(drawn);
    TileHeightCache cache;
    (void)builder.build(b, 17, &cache);  // fills the cache
    EXPECT_FALSE(cache.fine.empty());
    EXPECT_FALSE(cache.parent.empty());

    std::vector<std::pair<ChVector2i, double>> rut;
    for (int i = -30; i <= 30; ++i)
        for (int j = -3; j <= 3; ++j)
            rut.push_back({ChVector2i(i, j), -0.12});
    ruts->SetDeltas(rut);

    const Mesh cached = builder.build(b, 17, &cache);
    const Mesh full = builder.build(b, 17);
    EXPECT_EQ(cached.vertexData, full.vertexData);
    EXPECT_EQ(cached.minElevation, full.minElevation);
    EXPECT_EQ(cached.maxElevation, full.maxElevation);
}

TEST(ChPlanetQuadtree, SplitDistanceAndHorizonLimitTiles) {
    auto surface = moon::CreateSurface(15);
    const ChSiteFrame site = surface->MakeSiteFrame(30.75, 20.19);
    const ChVector3d cam = moon::Body().ToCartesian(30.75, 20.19, site.GetOriginElevation() + 2.0);
    auto count = [&](double split, double horizon) {
        ChPlanetQuadtree quadtree(surface, 2);
        quadtree.SetSplitDistance(split);
        quadtree.SetHorizonMargin(horizon);
        double t = 0;
        for (int k = 0; k < 10; ++k)
            quadtree.Update(cam, t += ChPlanetQuadtree::kTickS);
        int finest = 0;
        for (const ChTileMesh* m : quadtree.GetMeshes())
            finest = std::max(finest, m->level);
        EXPECT_EQ(finest, 17);  // full detail under the camera either way
        return quadtree.GetMeshes().size();
    };
    const size_t generous = count(1.79, -1.0);  // the previous behavior
    const size_t split_only = count(1.0, -1.0);
    const size_t both = count(1.0, 1000.0);
    EXPECT_LT(split_only, generous);
    EXPECT_LE(both, split_only);
    EXPECT_THROW(ChPlanetQuadtree(surface, 1).SetSplitDistance(0.0), std::invalid_argument);
}


TEST(CartesianTerrain, MeshCoordinatesSlopesAndSkirts) {
    auto surface = chrono_types::make_shared<ChPlanetSurface>(kBody);
    const ChSiteFrame site(kBody, 30, 20, 200);
    surface->SetSampler(chrono_types::make_shared<ChElevationFunction>([site](double lon, double lat) {
        const auto local = site.ToLocal(lon, lat, 200);
        return 200 + 2 * local.x() + 3 * local.y();
    }));
    const Cartesian::Boundary bounds{10, 20, 4, 4, site};
    const auto children = CoordinateTraits<Cartesian>::getChildBounds(bounds);
    EXPECT_EQ(children[0].centerX, 12);
    EXPECT_EQ(children[0].centerY, 22);
    EXPECT_EQ(children[2].centerX, 8);
    EXPECT_EQ(children[2].centerY, 18);
    EXPECT_DOUBLE_EQ(CoordinateTraits<Cartesian>::distanceToBounds(bounds, {10, 20, 110}, 60, 100, 0), 10);
    const TileMeshBuilder<Cartesian> builder(surface);
    const auto mesh = builder.build(bounds, 2);
    const int side = MeshTopology::gridSide(2);
    for (int row = 0; row < side; ++row) {
        for (int column = 0; column < side; ++column) {
            const size_t vertex = static_cast<size_t>(row) * side + column;
            const float* packed = &mesh.vertexData[vertex * ChTileMesh::kFloatsPerVertex];
            const double x = 6 + 8.0 * column / (side - 1), y = 16 + 8.0 * row / (side - 1);
            EXPECT_NEAR(mesh.centerX + packed[0], x, 1e-6);
            EXPECT_NEAR(mesh.centerY + packed[1], y, 1e-6);
            EXPECT_NEAR(mesh.centerZ + packed[2], 2 * x + 3 * y, 1e-6);
            EXPECT_NEAR(packed[3], -2, 1e-5);
            EXPECT_NEAR(packed[4], -3, 1e-5);
            EXPECT_NEAR(mesh.centerZ + packed[9], 2 * x + 3 * y, 1e-5);
        }
    }
    const size_t skirt = static_cast<size_t>(side) * side * ChTileMesh::kFloatsPerVertex;
    EXPECT_NEAR(mesh.vertexData[2] - mesh.vertexData[skirt + 2], builder.skirtDepthForLevel(2), 1e-4);
}

TEST(CartesianTerrain, CachedRebuildAndLod) {
    class Offset : public ChSurfaceFilter {
      public:
        double delta = 0;
        bool IsDynamic() const override { return true; }
        double Apply(double, double, double, double height) const override { return height + delta; }
    };
    auto surface = chrono_types::make_shared<ChPlanetSurface>(kBody);
    auto offset = chrono_types::make_shared<Offset>();
    surface->AddFilter(offset);
    const Cartesian::Boundary bounds{0, 0, 8, 8, ChSiteFrame(kBody, 179.9999, 20, 0)};
    TileHeightCache cache;
    const TileMeshBuilder<Cartesian> builder(surface);
    (void)builder.build(bounds, 2, &cache);
    offset->delta = -0.25;
    const auto cached = builder.build(bounds, 2, &cache);
    const auto full = builder.build(bounds, 2);
    EXPECT_EQ(cached.vertexData, full.vertexData);
    EXPECT_DOUBLE_EQ(cached.centerZ, full.centerZ);
    QuadtreeTile<Cartesian> tile(bounds, surface);
    const LodParams lod{16, 1000, kBody.GetRadius()};
    tile.updateLOD({0, 0, 10}, lod);
    EXPECT_EQ(tile.getMeshes().size(), 5u);
    offset->delta = -0.5;
    EXPECT_EQ(tile.invalidate(-8, -8, 8, 8, 10), 5);
    for (const auto& entry : tile.getMeshes())
        EXPECT_NEAR(entry.second.centerZ, -0.5, 1e-12);
    tile.updateLOD({1000, 0, 10}, lod);
    EXPECT_EQ(tile.getMeshes().size(), 1u);
    EXPECT_NEAR(tile.getElevation({0, 0}, 2), -0.5, 1e-12);
    // Splitting again builds all children before returning, with the current surface heights.
    tile.updateLOD({0, 0, 10}, lod);
    ASSERT_EQ(tile.getMeshes().size(), 5u);
    for (const auto& entry : tile.getMeshes()) {
        EXPECT_FALSE(entry.second.vertexData.empty());
        EXPECT_NEAR(entry.second.centerZ, -0.5, 1e-12);
    }
}

