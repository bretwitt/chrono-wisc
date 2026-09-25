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
// Tests of the meshing helpers: cutting holes in terrain meshes (ChSiteHoles)
// and meshing splatted spheres (ChSparseSdfGrid, as PlanetCRMWindow draws loose
// soil).
//
// =============================================================================

#include <cmath>
#include <map>
#include <tuple>
#include <vector>

#include "gtest/gtest.h"

#include "chrono/geometry/ChTriangleMeshConnected.h"

#include "chrono_planet/ChSiteFrame.h"
#include "chrono_planet/ChSiteHoles.h"
#include "chrono_planet/geometry/ChSdfMesher.h"

using namespace chrono;
using namespace chrono::planet;

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
// Sparse splat grid

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

