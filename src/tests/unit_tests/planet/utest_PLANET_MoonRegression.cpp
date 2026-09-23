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
// Regression test of the Moon preset. The reference values were captured from
// the module before its body, DEM and relief parameters became configurable,
// when the Moon was compiled in. The preset must keep reproducing them: physics
// grids and points, tile meshes, and quadtree heights,
// over synthetic DEMs (a global geographic base and a projected overlay with a
// nodata hole).
//
// =============================================================================

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "chrono_planet/lod/ChPlanetQuadtree.h"
#include "chrono_planet/lod/TileMeshBuilder.h"
#include "chrono_planet/planets/moon/ChMoon.h"

#include "planet_test_utils.h"

using namespace chrono::planet;
namespace golden = planet_test::golden;

namespace {

// Reference values: checksum of each grid, (bare, dem, dem3) per point, mesh checksums and geometry,
// and quadtree heights.
const std::map<std::string, std::vector<double>> kReference = {
    {"bare_grid site_fine", {643133.01367799821, 2571381.5461795563, 383.00264033075882, 382.18379573346186}},
    {"dem_grid site_fine", {3383739.0058649718, 13528915.045653475, 2013.5212052745533, 2012.0039952857164}},
    {"dem3_grid site_fine", {3383078.4572422663, 13526277.433870045, 2012.9818535409649, 2012.0957521364764}},
    {"bare_grid overlay_edge", {494892.65451211634, 1979494.0978258434, 296.0280879192909, 286.46625711861498}},
    {"dem_grid overlay_edge", {3247479.1220404506, 12984947.973491278, 1936.6896253003902, 1923.2111157474974}},
    {"dem3_grid overlay_edge", {3255311.3820499852, 13015486.061934356, 1949.2568504154065, 1923.6139039013519}},
    {"bare_grid overlay_hole", {396735.00380197744, 1584763.470043733, 370.20448761958409, 362.03370932989674}},
    {"dem_grid overlay_hole", {2191122.1709764278, 8752436.7574102208, 2009.3146971867823, 2018.495193126925}},
    {"dem3_grid overlay_hole", {2172144.4057885315, 8676609.0084218401, 1996.6463214625694, 1992.587669502233}},
    {"dem3_grid regional", {2950908.778297354, 11791368.195282871, 1765.9553113490476, 1208.2690844656383}},
    {"point site", {382.73644432525253, 2013.07882837223, 2012.5393100417521}},
    {"point east", {361.03039274592305, 2019.6004383579627, 1991.0628139340722}},
    {"point hole", {361.45769949019109, 2019.2481701642348, 1992.1229797758704}},
    {"point far", {127.21289554534593, 961.96707855626858, 961.70537573892022}},
    {"point west_wrap", {-24.458477032134343, 679.4739669752538, 679.42405340380117}},
    {"mesh_vertices root", {-261565.39924775343, 634280.18414511811, 140757.71875, -0.067727431654930115}},
    {"mesh_geom root", {1521552.1990630359, 677447.42444457987, 417967.09403792792, 342535.64079162944, 1360.9217495650519, 1893.3104994189925}},
    {"mesh_vertices l3", {-115241.11325192789, -421049.9222334458, 23934.640625, 0.0048130056820809841}},
    {"mesh_geom l3", {1391601.0552865611, 836157.0631011111, 623232.81196059601, 41864.649314084483, 1360.9217495650519, 2054.6583600509912}},
    {"mesh_vertices l6", {-40884.706107338308, -169388.27117718116, 2983.14697265625, 0.015991874039173126}},
    {"mesh_geom l6", {1401735.3504240338, 838090.80711696437, 598465.55730918224, 5276.8888698257979, 1501.4807398989797, 2122.6358961914666}},
    {"mesh_vertices l9", {-15226.37479064316, -61552.115375456779, 368.54318237304688, 0.13360701501369476}},
    {"mesh_geom l9", {1401015.2239537367, 837141.69495027058, 601596.52088649338, 670.56180047460782, 1992.6558508237358, 2020.7323339267168}},
    {"world zoom 0", {2011.8277555442392}},
    {"world zoom 5", {2020.4453668584313}},
    {"world zoom 10", {2020.4453668584313}},
    {"world zoom 15", {2020.435269799465}},
};

// Relative tolerance. Reordered sums differ in the last bits; anything larger is a real change.
constexpr double kRelTol = 1e-8;

void Expect(const std::string& key, const std::vector<double>& actual) {
    const auto it = kReference.find(key);
    ASSERT_NE(it, kReference.end()) << key;
    const auto& expected = it->second;
    ASSERT_EQ(expected.size(), actual.size()) << key;
    for (size_t i = 0; i < expected.size(); ++i)
        EXPECT_NEAR(actual[i], expected[i], kRelTol * std::max(1.0, std::abs(expected[i]))) << key << " [" << i << "]";
}

class MoonRegression : public ::testing::Test {
  protected:
    static void SetUpTestSuite() {
        base = planet_test::TempDir() + "/moon_base.tif";
        overlay = planet_test::TempDir() + "/moon_overlay.tif";
        golden::WriteDems(base, overlay, golden::kMoonRadius);
    }
    static std::vector<ChGeoTiffSource> Dems() { return {{base, 0, 4}, {overlay, 5, 30}}; }
    static std::string base, overlay;
};
std::string MoonRegression::base, MoonRegression::overlay;

}  // namespace

TEST_F(MoonRegression, PhysicsGrids) {
    auto bare = moon::CreateSurface(15);
    auto dem = moon::CreateSurface(Dems(), 15);
    auto dem3 = moon::CreateSurface(Dems(), 3);
    std::vector<double> h;
    for (const auto& g : golden::Grids()) {
        if (g.fine) {
            bare->GetElevationGrid(g.lon0, g.lat0, g.step, g.step, g.n, h);
            Expect(std::string("bare_grid ") + g.name, golden::Checksum(h));
            dem->GetElevationGrid(g.lon0, g.lat0, g.step, g.step, g.n, h);
            Expect(std::string("dem_grid ") + g.name, golden::Checksum(h));
        }
        dem3->GetElevationGrid(g.lon0, g.lat0, g.step, g.step, g.n, h);
        Expect(std::string("dem3_grid ") + g.name, golden::Checksum(h));
    }
}

TEST_F(MoonRegression, PhysicsPoints) {
    auto bare = moon::CreateSurface(15);
    auto dem = moon::CreateSurface(Dems(), 15);
    auto dem3 = moon::CreateSurface(Dems(), 3);
    for (const auto& p : golden::Points()) {
        Expect(std::string("point ") + p.name,
               {bare->GetElevation(p.lon, p.lat), dem->GetElevation(p.lon, p.lat), dem3->GetElevation(p.lon, p.lat)});
    }
}

TEST_F(MoonRegression, TileMeshes) {
    auto surface = moon::CreateSurface(Dems(), 15);
    for (const auto& t : golden::Tiles()) {
        const Mesh m = TileMeshBuilder<Spherical>(surface).build(Spherical::Boundary{t.lon, t.lat, t.half, t.half}, t.level);
        const std::string name = t.name;
        Expect("mesh_vertices " + name, golden::Checksum(m.vertexData));
        Expect("mesh_geom " + name,
               {m.centerX, m.centerY, m.centerZ, m.radius, m.minElevation, m.maxElevation});
    }
}

TEST_F(MoonRegression, QuadtreeElevation) {
    auto surface = moon::CreateSurface(Dems(), 15);
    ChPlanetQuadtree quadtree(surface, 1);
    const auto cam = surface->GetBody().ToCartesian(golden::kSiteLon, golden::kSiteLat, 3000.0);
    quadtree.Update(cam, 0.0);
    for (int zoom : {0, 5, 10, 15}) {
        const auto h = quadtree.GetElevation(golden::kSiteLon + 0.01, golden::kSiteLat - 0.02, zoom);
        ASSERT_TRUE(h.has_value());
        Expect("world zoom " + std::to_string(zoom), {*h});
    }
}
