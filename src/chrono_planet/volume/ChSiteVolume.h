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
// The ground of a work site as a signed distance field, so it can be dug and
// built up into shapes a height field cannot hold: trench walls, undercuts,
// spoil piles over cuts.
//
// =============================================================================

#ifndef CH_SITE_VOLUME_H
#define CH_SITE_VOLUME_H

#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <vector>

#include "chrono/core/ChFrame.h"
#include "chrono/core/ChVector3.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"

#include "chrono_planet/ChApiPlanet.h"
#include "chrono_planet/ChSiteFrame.h"
#include "chrono_planet/volume/ChSdfShape.h"

namespace chrono {
namespace planet {

class ChDeformationFilter;

/// @addtogroup planet_module
/// @{

/// The ground over a rectangle of a site frame as a signed distance field (SDF), negative in the soil and
/// positive in the space above it, which excavation edits as solids: Subtract digs a shape out, Add and Deposit
/// build ground up. Unlike the height field the rest of the terrain is, it holds any shape, such as a trench with
/// vertical or undercut walls, or a pile of spoil over a cut.
///
/// The field is sampled on a lattice of cubic voxels whose points sit at whole multiples of the voxel size in the
/// site frame, and stored only in a narrow band around the ground, in bricks of 8^3 points: bricks wholly in
/// the soil or the space above hold no samples. It starts as the terrain given by a height function, and the
/// region it covers is meant to be cut out of the drawn terrain and of the terrain physics collides with (see
/// ChPlanetVisualMesh::AddHole, ChSiteVolumeShapes). Edits stay within GetEditRegion(), a margin inside the
/// region, so where the volume meets the terrain around it both have the same ground.
///
/// Volumes changed by edits are measured on the lattice, each point standing for a voxel filled in proportion
/// to how deep in the soil it is, so what one edit removes and another puts back balance exactly. Meshers and
/// the top-surface publisher find what changed by versions (GetChangedBricks, PublishTopSurface).
///
/// Not thread-safe: edit and query it from one thread.
class CH_PLANET_API ChSiteVolume {
  public:
    /// Ground height (m) at site x/y. Called from the constructing thread only.
    using HeightFunction = std::function<double(double x, double y)>;

    /// Lattice points per brick along each axis.
    static constexpr int kBrick = 8;

    struct Params {
        double voxel = 0.02;   ///< voxel edge (m)
        double depth = 1.5;    ///< the volume reaches this far below the lowest ground of the region (m)
        double height = 1.5;   ///< and this far above the highest (m)
        double margin = 0.1;   ///< edits keep this far inside the region (m); at least two voxels are kept
    };

    /// Take over the ground of `region` from a height function.
    ChSiteVolume(const ChSiteRegion& region, const Params& params, HeightFunction ground);

    const Params& GetParams() const { return m_params; }
    double GetVoxelSize() const { return m_h; }
    /// The rectangle the volume covers, snapped out to whole voxels.
    const ChSiteRegion& GetRegion() const { return m_region; }
    /// The rectangle edits are confined to.
    const ChSiteRegion& GetEditRegion() const { return m_edit_region; }
    /// The rectangle to cut out of the terrain around the volume (ChSiteHoles, ChPlanetVisualMesh::AddHole): a voxel
    /// inside the region, so that the volume and the terrain both hold a narrow band of the same ground instead of
    /// leaving a crack, and holding the whole edit region.
    ChSiteRegion GetHole() const { return m_region.Inset(m_h); }
    /// Lowest and highest ground the volume can hold (m).
    double GetMinZ() const { return m_origin.z(); }
    double GetMaxZ() const { return m_origin.z() + (m_n[2] - 1) * m_h; }

    // -------------------------------------------------------------------------
    // Queries

    /// Signed distance (m) at a site point, trilinear between lattice points. Farther than a few voxels from the
    /// ground it is only a bound, and outside the volume's box it extends the nearest lattice column.
    double GetDistance(const ChVector3d& p) const;

    /// Unit normal of the ground nearest a point, from the gradient of the distance, pointing out of the soil.
    ChVector3d GetNormal(const ChVector3d& p) const;

    /// True if a point is in the soil.
    bool IsSolid(const ChVector3d& p) const { return GetDistance(p) < 0; }

    /// Height (m) of the topmost ground at site x/y, bilinear between lattice columns: what a vertical ray from
    /// above first meets. Outside the region, the nearest column's.
    double GetTopHeight(double x, double y) const;

    /// Height of the ground at site x/y before any edit (m).
    double GetInitialHeight(double x, double y) const;

    /// Volume of soil in the region (m^3), measured on the lattice, down to GetMinZ().
    double GetSoilVolume() const { return m_soil_volume; }

    // -------------------------------------------------------------------------
    // Edits. Each returns the volume of soil it changed (m^3, measured on the lattice) and confines itself to
    // GetEditRegion().

    /// Remove the soil inside a shape placed by `frame` (shape frame to site frame). Returns the volume removed.
    double Subtract(const ChSdfShape& shape, const ChFrame<>& frame);

    /// Remove the soil a shape sweeps through moving from one placement to another, the placements between
    /// interpolated no more than half a voxel apart. Returns the volume removed.
    double SubtractSwept(const ChSdfShape& shape, const ChFrame<>& from, const ChFrame<>& to);

    /// Fill a shape placed by `frame` with soil. Returns the volume added.
    double Add(const ChSdfShape& shape, const ChFrame<>& frame);

    /// Pour a volume of soil (m^3) onto the ground at site x/y as a cone at the angle of repose (rad), which
    /// also fills what lies under it, such as a trench. Returns the volume added, which differs from the one asked
    /// for by the lattice's resolution, or by more if the pile reaches the edge of the edit region.
    double Deposit(double x, double y, double volume, double repose_angle);

    /// Set the distance at the lattice points in a box (and the edit region) from a function of the site point,
    /// replacing what was there, for example with the ground a particle model computed there (vehicle::PlanetCRMTerrain).
    /// Returns the change of soil volume (m^3), positive if soil was added.
    double Assign(const ChAABB& box, const std::function<float(const ChVector3d& p)>& distance);

    // -------------------------------------------------------------------------
    // Changes

    /// Version of the latest edit that changed the field (0 before any, 1 for the initial ground).
    std::uint64_t GetVersion() const { return m_version; }

    /// Bricks changed after version `since`, or whose mesh changed: a brick meshes points of its neighbors, so
    /// they are included. A brick that came to hold no samples is included so its old mesh can be dropped.
    void GetChangedBricks(std::uint64_t since, std::vector<ChVector3i>& bricks) const;

    /// Bricks along each axis.
    ChVector3i GetNumBricks() const { return ChVector3i(m_nb[0], m_nb[1], m_nb[2]); }

    /// Number of bricks that hold samples.
    size_t GetNumAllocatedBricks() const { return m_pool.size() - m_free.size(); }

    /// Mesh the ground in one brick, with surface nets, in site coordinates: counter-clockwise triangles seen from
    /// the space above, vertex normals from the distance gradient. Faces on ground that edits moved by more than
    /// `disturbed` (m) get material index 1, the others 0. Meshes of neighboring bricks meet without gaps.
    /// Returns false, leaving the mesh empty, if the brick holds no ground.
    bool MeshBrick(const ChVector3i& brick, ChTriangleMeshConnected& mesh, double disturbed = 0.005) const;

    /// Write the changes to the top of the ground since the last call into a deformation filter, as the height
    /// change from the initial ground and the height after it (above the reference sphere of `site`, the frame
    /// the volume is in), so terrain drawn as a height field around the site, or from far away, shows what was
    /// dug and piled up. The filter's node spacing must equal the voxel size. Returns the nodes written.
    size_t PublishTopSurface(ChDeformationFilter& filter, const ChSiteFrame& site);

  private:
    using Brick = std::array<float, kBrick * kBrick * kBrick>;
    // Brick table entries for bricks without samples
    static constexpr std::int32_t kAbove = -1;
    static constexpr std::int32_t kBelow = -2;

    size_t BrickIndex(int bi, int bj, int bk) const { return size_t(bi) + size_t(m_nb[0]) * (size_t(bj) + size_t(m_nb[1]) * size_t(bk)); }
    size_t ColumnIndex(int i, int j) const { return size_t(i) + size_t(m_n[0]) * size_t(j); }
    ChVector3d PointPosition(int i, int j, int k) const { return m_origin + ChVector3d(i, j, k) * m_h; }

    // Sample at a lattice point, clamped to the lattice: above the top the space, below the bottom the soil
    float Sample(int i, int j, int k) const;
    // Initial distance at a lattice point from the initial ground
    float InitialSample(int i, int j, int k) const;
    // Soil fraction a sample stands for
    double Fill(float d) const;

    // Apply `op` (old sample, point) -> new sample over the lattice points in a box, within the edit region,
    // update the bricks, versions, soil volume and top heights, and return the change of soil volume.
    double Edit(const ChAABB& box, const std::function<float(float, const ChVector3d&)>& op);
    // The change of soil volume Edit would make, without making it
    double MeasureEdit(const ChAABB& box, const std::function<float(float, const ChVector3d&)>& op) const;
    // Lattice points in a box and the edit region: false if none
    bool EditRange(const ChAABB& box, int& i0, int& j0, int& k0, int& i1, int& j1, int& k1) const;
    // Recompute the top heights of lattice columns [i0, i1] x [j0, j1]
    void UpdateTops(int i0, int j0, int i1, int j1);
    // Top height of a lattice column
    double ColumnTop(int i, int j) const;

    Params m_params;
    double m_h;
    float m_band;               // samples are clamped to [-m_band, m_band]
    ChSiteRegion m_region;
    ChSiteRegion m_edit_region;
    ChVector3d m_origin;        // position of lattice point (0, 0, 0)
    int m_n[3];                 // lattice points along each axis
    int m_nb[3];                // bricks along each axis
    std::vector<float> m_ground;     // initial ground height of each lattice column
    std::vector<float> m_slope;      // 1 / sqrt(1 + |grad h|^2) of the initial ground, per column
    std::vector<float> m_top;        // current top height of each lattice column
    std::vector<std::int32_t> m_table;  // brick -> pool index, or kAbove / kBelow
    std::vector<std::uint64_t> m_brick_version;  // version of the last change to each brick's samples or mesh
    std::vector<Brick> m_pool;
    std::vector<std::int32_t> m_free;
    std::uint64_t m_version = 0;
    double m_soil_volume = 0;

    // Columns whose top changed since the last PublishTopSurface
    int m_dirty_i0, m_dirty_j0, m_dirty_i1, m_dirty_j1;
};

/// @} planet_module

}  // namespace planet
}  // namespace chrono

#endif
