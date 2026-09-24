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
// Tests of ChSiteVolume: the field it starts from, the volumes its edits
// change, where edits may reach, the meshes of its ground, and what it
// publishes to a deformation filter; of cutting holes in terrain meshes; of
// ChExcavationTool; and of ChSiteVolumeShapes.
//
// =============================================================================

#include <cmath>
#include <map>
#include <tuple>
#include <vector>

#include "gtest/gtest.h"

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/utils/ChConstants.h"

#include "chrono_planet/ChSiteFrame.h"
#include "chrono_planet/ChSiteHoles.h"
#include "chrono_planet/filters/ChDeformationFilter.h"
#include "chrono_planet/volume/ChSdfShape.h"
#include "chrono_planet/volume/ChExcavationTool.h"
#include "chrono_planet/volume/ChSiteVolume.h"
#include "chrono_planet/volume/ChSdfMesher.h"
#include "chrono_planet/volume/ChSiteVolumeShapes.h"
#include "chrono_planet/volume/ChSoilParticles.h"

using namespace chrono;
using namespace chrono::planet;

namespace {

const double kVoxel = 0.02;

ChSiteVolume::Params TestParams() {
    ChSiteVolume::Params params;
    params.voxel = kVoxel;
    params.depth = 1.0;
    params.height = 1.0;
    params.margin = 0.1;
    return params;
}

ChSiteVolume FlatVolume(double size = 3.0) {
    return ChSiteVolume(ChSiteRegion(-size / 2, -size / 2, size / 2, size / 2), TestParams(), [](double, double) { return 0.0; });
}

// Mesh every brick holding ground into one list of triangles
std::vector<std::array<ChVector3d, 3>> MeshAll(const ChSiteVolume& volume) {
    std::vector<std::array<ChVector3d, 3>> triangles;
    std::vector<ChVector3i> bricks;
    volume.GetChangedBricks(0, bricks);
    ChTriangleMeshConnected mesh;
    for (const auto& b : bricks) {
        if (!volume.MeshBrick(b, mesh))
            continue;
        const auto& v = mesh.GetCoordsVertices();
        for (const auto& f : mesh.GetIndicesVertices())
            triangles.push_back({v[f.x()], v[f.y()], v[f.z()]});
    }
    return triangles;
}

}  // namespace

TEST(ChSiteVolume, FlatGround) {
    ChSiteVolume volume = FlatVolume();
    EXPECT_NEAR(volume.GetTopHeight(0.3, -0.2), 0.0, 1e-6);
    EXPECT_NEAR(volume.GetDistance(ChVector3d(0.1, 0.1, 0.03)), 0.03, 1e-6);
    EXPECT_NEAR(volume.GetDistance(ChVector3d(0.1, 0.1, -0.025)), -0.025, 1e-6);
    EXPECT_TRUE(volume.IsSolid(ChVector3d(0, 0, -0.5)));
    EXPECT_FALSE(volume.IsSolid(ChVector3d(0, 0, 0.5)));
    EXPECT_NEAR(volume.GetNormal(ChVector3d(0.2, 0.2, 0.01)).z(), 1.0, 1e-6);
    // Only the bricks the ground passes through hold samples
    const ChVector3i nb = volume.GetNumBricks();
    EXPECT_LT(volume.GetNumAllocatedBricks(), size_t(nb.x()) * nb.y() * 2 + 1);
}

TEST(ChSiteVolume, SlopedGroundDistance) {
    const double slope = 0.3;
    ChSiteVolume volume(ChSiteRegion(-1, -1, 1, 1), TestParams(), [&](double x, double) { return slope * x; });
    const double x = 0.2, z = slope * x + 0.02;
    EXPECT_NEAR(volume.GetDistance(ChVector3d(x, 0.1, z)), 0.02 / std::sqrt(1 + slope * slope), 2e-4);
    EXPECT_NEAR(volume.GetTopHeight(x, 0.1), slope * x, 1e-3);
}

TEST(ChSiteVolume, SubtractSphereVolume) {
    ChSiteVolume volume = FlatVolume();
    const double before = volume.GetSoilVolume();
    const double r = 0.4;
    const double removed = volume.Subtract(ChSdfSphere(r), ChFrame<>(ChVector3d(0, 0, 0)));
    const double half_sphere = 2.0 / 3.0 * CH_PI * r * r * r;
    EXPECT_NEAR(removed, half_sphere, 0.02 * half_sphere);
    EXPECT_NEAR(volume.GetSoilVolume(), before - removed, 1e-9);
    EXPECT_NEAR(volume.GetTopHeight(0, 0), -r, 0.01);
    EXPECT_NEAR(volume.GetTopHeight(0.5, 0), 0.0, 1e-6);
    // Removing it again removes nothing
    EXPECT_NEAR(volume.Subtract(ChSdfSphere(r), ChFrame<>(ChVector3d(0, 0, 0))), 0.0, 1e-9);
}

TEST(ChSiteVolume, SweptCapsuleVolume) {
    // A ball dragged through the soil well below the ground cuts a tunnel: a cylinder and two half balls
    ChSiteVolume volume = FlatVolume();
    const double r = 0.15, length = 0.8, depth = -0.5;
    const double removed = volume.SubtractSwept(ChSdfSphere(r), ChFrame<>(ChVector3d(-length / 2, 0, depth)), ChFrame<>(ChVector3d(length / 2, 0, depth)));
    const double expected = CH_PI * r * r * length + 4.0 / 3.0 * CH_PI * r * r * r;
    EXPECT_NEAR(removed, expected, 0.03 * expected);
    // The ground above the tunnel is untouched: an undercut a height field could not hold
    EXPECT_NEAR(volume.GetTopHeight(0, 0), 0.0, 1e-6);
    EXPECT_FALSE(volume.IsSolid(ChVector3d(0, 0, depth)));
    EXPECT_TRUE(volume.IsSolid(ChVector3d(0, 0, depth + r + 0.05)));
}

TEST(ChSiteVolume, DepositConservesVolume) {
    ChSiteVolume volume = FlatVolume(4.0);
    const double initial = volume.GetSoilVolume();
    const double removed = volume.Subtract(ChSdfBox(ChVector3d(0.6, 0.4, 0.6)), ChFrame<>(ChVector3d(-0.8, 0, 0)));
    EXPECT_NEAR(removed, 0.6 * 0.4 * 0.3, 0.03 * 0.072);
    const double added = volume.Deposit(0.8, 0, removed, 35 * CH_DEG_TO_RAD);
    EXPECT_NEAR(added, removed, 0.05 * removed);
    EXPECT_NEAR(volume.GetSoilVolume(), initial - removed + added, 1e-9);

    // On flat ground the pile is a cone of the angle of repose
    const double slope = std::tan(35 * CH_DEG_TO_RAD);
    const double height = std::cbrt(3 * removed * slope * slope / CH_PI);
    EXPECT_NEAR(volume.GetTopHeight(0.8, 0), height, 0.1 * height);
}

TEST(ChSiteVolume, DepositFillsTrench) {
    ChSiteVolume volume = FlatVolume();
    const double removed = volume.Subtract(ChSdfBox(ChVector3d(0.4, 1.0, 0.4)), ChFrame<>(ChVector3d(0, 0, 0)));
    // Pouring the soil back into the trench fills it before piling up
    volume.Deposit(0, 0, removed, 35 * CH_DEG_TO_RAD);
    EXPECT_GT(volume.GetTopHeight(0, 0), -0.05);
}

TEST(ChSiteVolume, EditsStayInEditRegion) {
    ChSiteVolume volume = FlatVolume(2.0);
    const ChSiteRegion& edit = volume.GetEditRegion();
    EXPECT_GT(edit.min_x, volume.GetRegion().min_x);
    // A sphere over the corner of the region digs only inside the edit region
    volume.Subtract(ChSdfSphere(0.5), ChFrame<>(ChVector3d(volume.GetRegion().max_x, volume.GetRegion().max_y, 0)));
    const double outside_x = 0.5 * (edit.max_x + volume.GetRegion().max_x);
    EXPECT_NEAR(volume.GetTopHeight(outside_x, outside_x), 0.0, 1e-6);
    EXPECT_LT(volume.GetTopHeight(edit.max_x - 0.02, edit.max_y - 0.02), -0.1);
}

TEST(ChSiteVolume, MeshIsClosedAndFacesOut) {
    // Edits placed off the lattice: where the ground passes exactly through a lattice point, the cells around it
    // put their vertices on that point, and welding by position below would join them
    ChSiteVolume volume = FlatVolume(2.0);
    volume.Subtract(ChSdfSphere(0.3517), ChFrame<>(ChVector3d(0.1013, -0.0471, 0.0537)));
    volume.SubtractSwept(ChSdfSphere(0.1031), ChFrame<>(ChVector3d(-0.4043, 0.3017, -0.3029)), ChFrame<>(ChVector3d(0.2957, 0.3083, -0.3511)));
    volume.Deposit(-0.4011, -0.3987, 0.02, 35 * CH_DEG_TO_RAD);
    const auto triangles = MeshAll(volume);
    ASSERT_GT(triangles.size(), 1000u);

    // Weld vertices by position; bricks mesh the cells on their borders alike, so seams weld exactly
    std::map<std::tuple<long, long, long>, int> ids;
    auto id = [&](const ChVector3d& p) {
        const auto key = std::make_tuple(std::lround(p.x() * 1e6), std::lround(p.y() * 1e6), std::lround(p.z() * 1e6));
        return ids.emplace(key, static_cast<int>(ids.size())).first->second;
    };
    std::map<std::pair<int, int>, int> directed;
    const ChSiteRegion inner = volume.GetRegion().Inset(3 * kVoxel);
    for (const auto& t : triangles) {
        int v[3] = {id(t[0]), id(t[1]), id(t[2])};
        for (int e = 0; e < 3; ++e)
            ++directed[{v[e], v[(e + 1) % 3]}];
        // Faces point out of the soil
        const ChVector3d n = Vcross(t[1] - t[0], t[2] - t[0]);
        const ChVector3d c = (t[0] + t[1] + t[2]) / 3;
        if (n.Length() > 1e-12) {
            EXPECT_GT(volume.GetDistance(c + n.GetNormalized() * (0.5 * kVoxel)), volume.GetDistance(c - n.GetNormalized() * (0.5 * kVoxel)));
        }
    }
    // Away from the region's edges, every edge is shared by exactly two faces, in opposite directions
    std::map<int, ChVector3d> position;
    for (const auto& t : triangles)
        for (int k = 0; k < 3; ++k)
            position[id(t[k])] = t[k];
    int open = 0;
    for (const auto& [edge, count] : directed) {
        const ChVector3d& a = position[edge.first];
        const ChVector3d& b = position[edge.second];
        if (edge.first == edge.second)
            continue;  // a sliver face
        if (!inner.Contains(a.x(), a.y()) || !inner.Contains(b.x(), b.y()))
            continue;
        const auto back = directed.find({edge.second, edge.first});
        if (count != 1 || back == directed.end() || back->second != 1)
            ++open;
    }
    EXPECT_EQ(open, 0);
}

TEST(ChSiteVolume, MeshMarksDisturbedGround) {
    ChSiteVolume volume = FlatVolume(2.0);
    volume.Subtract(ChSdfSphere(0.3), ChFrame<>(ChVector3d(0, 0, 0)));
    std::vector<ChVector3i> bricks;
    volume.GetChangedBricks(0, bricks);
    ChTriangleMeshConnected mesh;
    int disturbed = 0, undisturbed_far = 0, disturbed_far = 0;
    for (const auto& b : bricks) {
        if (!volume.MeshBrick(b, mesh))
            continue;
        const auto& v = mesh.GetCoordsVertices();
        const auto& f = mesh.GetIndicesVertices();
        const auto& m = mesh.GetIndicesMaterials();
        ASSERT_EQ(m.size(), f.size());
        for (size_t n = 0; n < f.size(); ++n) {
            const ChVector3d c = (v[f[n].x()] + v[f[n].y()] + v[f[n].z()]) / 3;
            const double r = std::hypot(c.x(), c.y());
            if (r < 0.2)
                disturbed += m[n] == 1;
            if (r > 0.6) {
                undisturbed_far += m[n] == 0;
                disturbed_far += m[n] == 1;
            }
        }
    }
    EXPECT_GT(disturbed, 0);
    EXPECT_GT(undisturbed_far, 0);
    EXPECT_EQ(disturbed_far, 0);
}

TEST(ChSiteVolume, ChangedBricks) {
    ChSiteVolume volume = FlatVolume(3.0);
    std::vector<ChVector3i> all, changed;
    volume.GetChangedBricks(0, all);
    EXPECT_EQ(all.size(), volume.GetNumAllocatedBricks());
    const std::uint64_t version = volume.GetVersion();
    volume.Subtract(ChSdfSphere(0.1), ChFrame<>(ChVector3d(0.5, 0.5, 0)));
    EXPECT_GT(volume.GetVersion(), version);
    volume.GetChangedBricks(version, changed);
    EXPECT_FALSE(changed.empty());
    EXPECT_LT(changed.size(), all.size() / 4);
    // An edit that changes nothing leaves the version alone
    const std::uint64_t after = volume.GetVersion();
    volume.Subtract(ChSdfSphere(0.1), ChFrame<>(ChVector3d(0.5, 0.5, 0.5)));
    EXPECT_EQ(volume.GetVersion(), after);
}

TEST(ChSiteVolume, PublishTopSurface) {
    const ChSiteFrame site(1737400.0, 10.0, 20.0, 100.0);
    ChSiteVolume volume = FlatVolume(2.0);
    auto filter = std::make_shared<ChDeformationFilter>(site, kVoxel);
    EXPECT_EQ(volume.PublishTopSurface(*filter, site), 0u);
    volume.Subtract(ChSdfBox(ChVector3d(0.4, 0.4, 0.4)), ChFrame<>(ChVector3d(0, 0, 0)));
    EXPECT_GT(volume.PublishTopSurface(*filter, site), 0u);
    ChGeoRegion region;
    EXPECT_GT(filter->GetChanges(0, region), 0u);
    double lon, lat;
    site.ToLonLat(0, 0, lon, lat);
    EXPECT_NEAR(filter->GetDelta(lon, lat, 0.0), -0.2, 0.01);
    // A filter of another spacing is refused
    ChDeformationFilter coarse(site, 0.05);
    EXPECT_THROW(volume.PublishTopSurface(coarse, site), std::invalid_argument);
}

// -----------------------------------------------------------------------------
// Holes cut into terrain meshes

namespace {
// A flat grid of n x n unit squares over [0, n]^2, two faces each, with normals, UVs and a material per face
ChTriangleMeshConnected Grid(int n) {
    ChTriangleMeshConnected mesh;
    auto& v = mesh.GetCoordsVertices();
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i) {
            v.push_back(ChVector3d(i, j, 0.1 * i));
            mesh.GetCoordsNormals().push_back(ChVector3d(0, 0, 1));
            mesh.GetCoordsUV().push_back(ChVector2d(i, j));
        }
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            const int a = j * (n + 1) + i, b = a + 1, c = a + n + 1, d = c + 1;
            for (const ChVector3i& f : {ChVector3i(a, b, d), ChVector3i(a, d, c)}) {
                mesh.GetIndicesVertices().push_back(f);
                mesh.GetIndicesNormals().push_back(f);
                mesh.GetIndicesUV().push_back(f);
                mesh.GetIndicesMaterials().push_back((i + j) % 2);
            }
        }
    return mesh;
}

double ProjectedArea(const ChTriangleMeshConnected& mesh) {
    double area = 0;
    const auto& v = mesh.GetCoordsVertices();
    for (const auto& f : mesh.GetIndicesVertices()) {
        const ChVector3d a = v[f[1]] - v[f[0]], b = v[f[2]] - v[f[0]];
        area += 0.5 * std::abs(a.x() * b.y() - a.y() * b.x());
    }
    return area;
}
}  // namespace

TEST(ChSiteHoles, CutsExactlyTheRectangle) {
    ChTriangleMeshConnected mesh = Grid(10);
    const ChSiteRegion hole(2.5, 3.25, 7.5, 6.75);
    EXPECT_GT(CutSiteHoles(mesh, {hole}), 0u);
    EXPECT_NEAR(ProjectedArea(mesh), 100.0 - 5.0 * 3.5, 1e-9);
    const auto& v = mesh.GetCoordsVertices();
    const size_t nf = mesh.GetIndicesVertices().size();
    ASSERT_EQ(mesh.GetIndicesNormals().size(), nf);
    ASSERT_EQ(mesh.GetIndicesUV().size(), nf);
    ASSERT_EQ(mesh.GetIndicesMaterials().size(), nf);
    for (size_t f = 0; f < nf; ++f) {
        const ChVector3i& face = mesh.GetIndicesVertices()[f];
        const ChVector3d c = (v[face[0]] + v[face[1]] + v[face[2]]) / 3;
        EXPECT_FALSE(hole.Inset(1e-9).Contains(c.x(), c.y()));
        // New vertices stay on the face's plane and carry its interpolated UVs
        for (int k = 0; k < 3; ++k) {
            EXPECT_NEAR(v[face[k]].z(), 0.1 * v[face[k]].x(), 1e-9);
            const ChVector2d& uv = mesh.GetCoordsUV()[mesh.GetIndicesUV()[f][k]];
            EXPECT_NEAR(uv.x(), v[face[k]].x(), 1e-9);
            EXPECT_NEAR(uv.y(), v[face[k]].y(), 1e-9);
        }
    }
    // A mesh clear of the hole is untouched
    ChTriangleMeshConnected other = Grid(2);
    EXPECT_EQ(CutSiteHoles(other, {ChSiteRegion(5, 5, 6, 6)}), 0u);
    EXPECT_EQ(other.GetIndicesVertices().size(), 8u);
}

// -----------------------------------------------------------------------------
// Excavation tool

TEST(ChExcavationTool, CuttingForceScaling) {
    ChExcavationTool::Soil soil;
    soil.cohesion = 0;
    const double f1 = ChExcavationTool::CuttingForce(soil, 0.5, 0.1, 45 * CH_DEG_TO_RAD);
    const double f2 = ChExcavationTool::CuttingForce(soil, 0.5, 0.2, 45 * CH_DEG_TO_RAD);
    const double f3 = ChExcavationTool::CuttingForce(soil, 1.0, 0.1, 45 * CH_DEG_TO_RAD);
    EXPECT_GT(f1, 0);
    EXPECT_NEAR(f2 / f1, 4.0, 1e-9);  // gamma term: d^2
    EXPECT_NEAR(f3 / f1, 2.0, 1e-9);  // width
    soil.cohesion = 1000;
    soil.bulk_density = 0;
    const double c1 = ChExcavationTool::CuttingForce(soil, 0.5, 0.1, 45 * CH_DEG_TO_RAD);
    const double c2 = ChExcavationTool::CuttingForce(soil, 0.5, 0.2, 45 * CH_DEG_TO_RAD);
    EXPECT_NEAR(c2 / c1, 2.0, 1e-9);  // cohesion term: d
    EXPECT_EQ(ChExcavationTool::CuttingForce(soil, 0.5, 0.0, 45 * CH_DEG_TO_RAD), 0.0);
}

TEST(ChExcavationTool, DigCarryDump) {
    auto volume = std::make_shared<ChSiteVolume>(FlatVolume(4.0));
    const double initial = volume->GetSoilVolume();
    ChSystemNSC sys;
    auto body = chrono_types::make_shared<ChBody>();
    body->SetFixed(true);
    sys.AddBody(body);
    // A box scoop 0.3 m wide, 0.2 m long and 0.2 m tall, dragged 1 m through the ground with its bottom 0.1 m deep
    auto cutter = std::make_shared<ChSdfBox>(ChVector3d(0.2, 0.3, 0.2));
    ChExcavationTool tool(volume, body, cutter);
    const double z = -0.1 + 0.1;
    body->SetPos(ChVector3d(-0.6, 0, z));
    tool.Update();
    const int steps = 200;
    for (int s = 1; s <= steps; ++s) {
        body->SetPos(ChVector3d(-0.6 + 1.0 * s / steps, 0, z));
        body->SetPosDt(ChVector3d(0.5, 0, 0));
        tool.Update();
    }
    // Cut: the trench under the ground, 1.2 m long counting the scoop's own length
    const double expected = 0.3 * 0.1 * 1.2;
    EXPECT_NEAR(tool.GetExcavatedVolume(), expected, 0.05 * expected);
    EXPECT_NEAR(tool.GetPayloadVolume(), tool.GetExcavatedVolume(), 1e-12);
    EXPECT_NEAR(volume->GetSoilVolume(), initial - tool.GetExcavatedVolume(), 1e-9);
    // The soil resists the motion and presses the scoop down, plus the payload's weight
    EXPECT_NEAR(tool.GetCuttingDepth(), 0.1, 0.02);
    EXPECT_LT(tool.GetResistanceForce().x(), 0);
    EXPECT_LT(tool.GetResistanceForce().z(), 0);

    // Dumped elsewhere, the soil comes back
    const double carried = tool.GetPayloadVolume();
    EXPECT_NEAR(tool.Dump(0.0, 1.2), carried, 1e-12);
    EXPECT_EQ(tool.GetPayloadVolume(), 0.0);
    EXPECT_NEAR(volume->GetSoilVolume(), initial, 0.05 * carried);
    EXPECT_GT(volume->GetTopHeight(0.0, 1.2), 0.05);
}

TEST(ChExcavationTool, FullToolPushesSoilAhead) {
    auto volume = std::make_shared<ChSiteVolume>(FlatVolume(4.0));
    const double initial = volume->GetSoilVolume();
    ChSystemNSC sys;
    auto body = chrono_types::make_shared<ChBody>();
    body->SetFixed(true);
    sys.AddBody(body);
    ChExcavationTool tool(volume, body, std::make_shared<ChSdfBox>(ChVector3d(0.2, 0.3, 0.2)));
    tool.SetCapacity(0.01);
    tool.EnableForces(false);
    body->SetPos(ChVector3d(-0.8, 0, 0));
    tool.Update();
    for (int s = 1; s <= 200; ++s) {
        body->SetPos(ChVector3d(-0.8 + 1.0 * s / 200, 0, 0));
        body->SetPosDt(ChVector3d(0.5, 0, 0));
        tool.Update();
    }
    EXPECT_NEAR(tool.GetPayloadVolume(), 0.01, 1e-12);
    EXPECT_GT(tool.GetSpilledVolume(), 0.0);
    // Spilled soil stays on the site: the volume lost is what the tool carries, less what still waits to be piled
    EXPECT_NEAR(volume->GetSoilVolume(), initial - tool.GetPayloadVolume(), 0.1 * tool.GetExcavatedVolume());
}

// -----------------------------------------------------------------------------
// Visual and collision shapes

TEST(ChSiteVolumeShapes, RebuildsChangedBricks) {
    ChSystemNSC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    auto volume = std::make_shared<ChSiteVolume>(FlatVolume(2.0));
    ChSiteVolumeShapes shapes(&sys, volume);
    shapes.SetDisturbedMaterial(chrono_types::make_shared<ChVisualMaterial>());
    shapes.EnableCollision(chrono_types::make_shared<ChContactMaterialNSC>(), 14);
    const size_t bodies = sys.GetBodies().size();
    const size_t first = shapes.Update();
    EXPECT_EQ(first, volume->GetNumAllocatedBricks());
    EXPECT_GT(shapes.GetNumTriangles(), 0u);
    EXPECT_EQ(shapes.GetBody()->GetVisualModel()->GetNumShapes(), shapes.GetNumShapes());
    EXPECT_EQ(sys.GetBodies().size(), bodies + shapes.GetNumShapes());
    EXPECT_EQ(shapes.Update(), 0u);

    volume->Subtract(ChSdfSphere(0.2), ChFrame<>(ChVector3d(0.1, 0.1, 0)));
    const size_t rebuilt = shapes.Update();
    EXPECT_GT(rebuilt, 0u);
    EXPECT_LT(rebuilt, first / 2);
    EXPECT_EQ(sys.GetBodies().size(), bodies + shapes.GetNumShapes());
    for (const auto& body : sys.GetBodies())
        if (body != shapes.GetBody() && body->IsCollisionEnabled())
            EXPECT_EQ(body->GetCollisionModel()->GetFamily(), 14);
}

// -----------------------------------------------------------------------------
// Sparse splat grid and soil particles

namespace {
// Count the edges, welded by position, that faces do not use as often in one direction as in the other: a closed
// surface has none. (Where two blobs touch through a neck one cell wide, surface nets puts four faces on an edge,
// two each way: closed, though not manifold.)
int OpenEdges(const ChTriangleMeshConnected& mesh) {
    std::map<std::tuple<long, long, long>, int> ids;
    auto id = [&](const ChVector3d& p) {
        const auto key = std::make_tuple(std::lround(p.x() * 1e6), std::lround(p.y() * 1e6), std::lround(p.z() * 1e6));
        return ids.emplace(key, static_cast<int>(ids.size())).first->second;
    };
    std::map<std::pair<int, int>, int> directed;
    const auto& v = mesh.GetCoordsVertices();
    for (const auto& f : mesh.GetIndicesVertices()) {
        const int a[3] = {id(v[f[0]]), id(v[f[1]]), id(v[f[2]])};
        for (int e = 0; e < 3; ++e)
            if (a[e] != a[(e + 1) % 3])
                ++directed[{a[e], a[(e + 1) % 3]}];
    }
    int open = 0;
    for (const auto& [edge, count] : directed) {
        const auto back = directed.find({edge.second, edge.first});
        if (back == directed.end() || back->second != count)
            ++open;
    }
    return open;
}
}  // namespace

TEST(ChSparseSdfGrid, ClosedAcrossBricks) {
    ChSparseSdfGrid grid(0.02, 0.1f);
    grid.Begin();
    // A clump of spheres spanning several bricks, placed off the lattice
    for (int n = 0; n < 60; ++n)
        grid.SplatSphere(ChVector3d(0.0137 * n, 0.3 * std::sin(0.2 * n) + 0.0071, 0.05 * std::cos(0.3 * n) - 0.0033), 0.025, 0.02);
    EXPECT_GT(grid.Remesh(), 0u);
    ChTriangleMeshConnected mesh;
    grid.AppendMeshes(mesh);
    ASSERT_GT(mesh.GetNumTriangles(), 100u);
    EXPECT_GT(grid.GetNumBricks(), 4u);
    EXPECT_EQ(OpenEdges(mesh), 0);

    // The same splats again change nothing to mesh
    grid.Begin();
    for (int n = 0; n < 60; ++n)
        grid.SplatSphere(ChVector3d(0.0137 * n, 0.3 * std::sin(0.2 * n) + 0.0071, 0.05 * std::cos(0.3 * n) - 0.0033), 0.025, 0.02);
    EXPECT_EQ(grid.Remesh(), 0u);
    // Nothing splatted: the bricks go
    grid.Begin();
    grid.Remesh();
    EXPECT_EQ(grid.GetNumBricks(), 0u);
}

TEST(ChSoilParticles, CarryPourLand) {
    auto volume = std::make_shared<ChSiteVolume>(FlatVolume(4.0));
    const double initial = volume->GetSoilVolume();
    ChSystemNSC sys;
    auto body = chrono_types::make_shared<ChBody>();
    body->SetFixed(true);
    sys.AddBody(body);
    const ChVector3d size(0.2, 0.3, 0.2);
    ChExcavationTool tool(volume, body, std::make_shared<ChSdfBox>(size));
    auto particles = std::make_shared<ChSoilParticles>(&sys, volume);
    tool.SetParticles(particles, ChAABB(-0.5 * size, 0.5 * size));
    tool.EnableForces(false);
    const double dt = 2e-3;
    auto step = [&]() {
        sys.SetChTime(sys.GetChTime() + dt);
        tool.Update();
        particles->Advance(dt);
    };

    // Drag the scoop through the ground: what it cuts shows as particles in it
    body->SetPos(ChVector3d(-0.6, 0, 0));
    tool.Update();
    for (int s = 1; s <= 200; ++s) {
        body->SetPos(ChVector3d(-0.6 + 0.6 * s / 200, 0, 0));
        body->SetPosDt(ChVector3d(0.3, 0, 0));
        step();
    }
    const double v = particles->GetParticleVolume();
    EXPECT_GT(tool.GetPayloadVolume(), 10 * v);
    EXPECT_EQ(particles->GetNumCarried(), static_cast<size_t>(std::floor(tool.GetPayloadVolume() / v + 1e-9)));
    particles->UpdateMesh();
    EXPECT_GT(particles->GetMesh()->GetNumTriangles(), 0u);

    // Lift it over fresh ground and tip it: the load pours out and falls, still in view
    body->SetPosDt(VNULL);
    body->SetPos(ChVector3d(1.0, 0.8, 0.8));
    step();
    EXPECT_EQ(particles->GetNumFalling(), 0u);
    const double load = tool.GetPayloadVolume();
    // Tipped over in half a second and held there: a load takes about the pour time to empty past vertical
    bool seen_falling = false;
    for (int s = 0; s < 1000; ++s) {
        body->SetRot(QuatFromAngleY(-std::min(2.0, 2.0 * s * dt / 0.5)));
        step();
        seen_falling = seen_falling || particles->GetNumFalling() > 0;
        if (s % 50 == 0) {
            particles->UpdateMesh();
            EXPECT_GT(particles->GetMesh()->GetNumTriangles(), 0u);  // in the tool, pouring, or falling
        }
    }
    EXPECT_TRUE(seen_falling);
    EXPECT_LT(tool.GetPayloadVolume(), v);
    EXPECT_NEAR(tool.GetDumpedVolume(), load - tool.GetPayloadVolume(), 1e-12);
    // Soil is conserved at every stage: in the ground, in the tool, or in flight
    EXPECT_NEAR(volume->GetSoilVolume() + tool.GetPayloadVolume() + particles->GetFreeVolume() + particles->GetLostVolume(), initial, 1e-6);

    // It lands and becomes ground where it fell
    for (int s = 0; s < 1500 && particles->GetFreeVolume() > 0; ++s)
        step();
    particles->Merge();
    EXPECT_EQ(particles->GetFreeVolume(), 0.0);
    // Piles hold their volume to the precision their apex is found to
    EXPECT_NEAR(particles->GetLostVolume(), 0.0, 1e-4 * load);
    EXPECT_NEAR(volume->GetSoilVolume() + tool.GetPayloadVolume() + particles->GetLostVolume(), initial, 1e-6);
    // A pile below the tool's lip, which the tipped tool poured over
    double highest = 0;
    for (double x = 0.5; x <= 1.5; x += 0.05)
        for (double y = 0.3; y <= 1.3; y += 0.05)
            highest = std::max(highest, volume->GetTopHeight(x, y));
    EXPECT_GT(highest, 0.05);
}
