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
// Regolith lofted into ballistic flight on an airless body, for example by
// rover wheels, and its optical density on a voxel grid for rendering.
//
// =============================================================================

#ifndef CH_DUST_FIELD_H
#define CH_DUST_FIELD_H

#include <climits>
#include <cstdint>
#include <functional>
#include <random>
#include <string>
#include <vector>

#include "chrono/core/ChVector3.h"

#include "chrono_planet/ChApiPlanet.h"

namespace chrono {
namespace planet {

/// @addtogroup planet_module
/// @{

/// Regolith grains in ballistic flight over the ground of an airless body, in a z-up Cartesian frame such as a
/// ChSiteFrame. With no gas to drag them, each grain follows a parabola from its release until it lands, so the
/// field holds only release states and evaluates positions in closed form at any time: emission can follow the
/// physics step while positions are needed only when a sensor looks.
///
/// Grains are tracked as super-particles, each standing for a mass of grains of one size bin. The optical
/// density of the field is binned on a voxel grid (UpdateGrid) as an extinction coefficient: a bin's grains
/// block 3 Q / (4 rho r) square meters of light per kilogram, so the fine grains, a small share of the mass,
/// set how thick the dust looks. The grid also carries, per voxel, the transmittance of the dust toward the Sun
/// and the Sun's visibility past the terrain, so a renderer can light the dust with a single lookup.
///
/// References: Hersh et al. 2012 (Am. J. Phys. 80, 452) for the ballistic dust trails of the Apollo Lunar Roving
/// Vehicle; Hapke 2012 (Theory of Reflectance and Emittance Spectroscopy) for leaving diffraction in the direct
/// beam, which makes Q about 1 for grains much larger than the wavelength.
class CH_PLANET_API ChDustField {
  public:
    /// Ground height (m) at a point (x, y) of the frame. Called from the calling thread only.
    using HeightFunction = std::function<double(double x, double y)>;

    /// A grain size bin: grain radius (m) and its share of the lofted mass.
    struct SizeBin {
        double radius;
        double mass_fraction;
    };

    /// Grain and soil properties.
    struct Params {
        double gravity = 1.62;               ///< gravitational acceleration (m/s^2), along -z
        double grain_density = 3100.0;       ///< density of a grain (kg/m^3)
        double extinction_efficiency = 1.0;  ///< Q, extinction cross section over geometric cross section
        /// Grain radii (m) and mass shares. The default is a coarse split of the lunar soil size distribution:
        /// about a sixth of the mass finer than 20 um across and roughly half coarser than 100 um.
        std::vector<SizeBin> bins = {{5e-6, 0.15}, {25e-6, 0.40}, {100e-6, 0.45}};
        double max_flight_time = 20.0;  ///< a grain still aloft after this long is dropped (s)
        size_t max_particles = 400000;  ///< emission stops while this many super-particles are aloft
        std::uint64_t seed = 1;         ///< random number seed
    };

    /// Loosened soil a wheel throws up.
    struct WheelEmission {
        double bulk_density = 1500.0;  ///< density of the loose surface soil (kg/m^3)
        /// Depth of loose soil the tread carries away per unit of rim travel (m). The lofted mass rate is
        /// bulk_density * width * loose_depth * (rim speed + slip_gain * slip speed).
        double loose_depth = 2e-4;
        double slip_gain = 4.0;  ///< extra loosening per unit of slip speed, relative to rim speed
        /// The rim releases soil from where it leaves the ground up to this angle from the bottom of the wheel
        /// (rad); 90 deg is level with the hub behind it. A fender lowers it.
        double max_release_angle = 1.5708;
        double min_speed_fraction = 0.3;  ///< released grains move at a random fraction of the rim velocity,
        double max_speed_fraction = 1.0;  ///< between these
        double spread = 0.1;              ///< random velocity added, relative to the rim speed (1 sigma per axis)
        double contact_tolerance = 0.02;  ///< the wheel emits while its bottom is within this of the ground (m)
        /// Super-particles per wheel and second, over all size bins. Each bin gets a share in proportion to the light
        /// its grains block, so each super-particle weighs about the same in the rendered dust.
        double particles_per_second = 20000.0;
    };

    /// State of one wheel.
    struct WheelState {
        ChVector3d center;    ///< wheel center (m)
        ChVector3d axle;      ///< unit spin axis
        ChVector3d velocity;  ///< velocity of the center (m/s)
        double omega;         ///< spin rate about the axle (rad/s, right-handed)
        double radius;        ///< wheel radius (m)
        double width;         ///< wheel width (m)
    };

    /// A super-particle: its release state and the mass of grains it stands for.
    struct Particle {
        ChVector3d pos0;  ///< position at release (m)
        ChVector3d vel0;  ///< velocity at release (m/s)
        double time0;     ///< release time (s)
        float mass;       ///< mass (kg)
        int bin;          ///< size bin
    };

    /// Size and resolution of the voxel grid.
    struct GridSpec {
        double voxel = 0.08;  ///< voxel edge (m)
        int nx = 160;         ///< voxels along x
        int ny = 160;         ///< voxels along y
        int nz = 40;          ///< voxels along z
        double below = 0.5;   ///< the grid starts this far below the center given to UpdateGrid (m)
        /// Beyond the grid, the Sun's visibility past the terrain is traced this far (m), so ground outside the
        /// grid still shades the dust.
        double shadow_distance = 40.0;
    };

    /// The binned field. Voxel (i, j, k) is centered at origin + voxel * (i + 1/2, j + 1/2, k + 1/2) and stored
    /// at index i + nx * (j + ny * k).
    struct Grid {
        ChVector3d origin;
        double voxel = 0;
        int nx = 0, ny = 0, nz = 0;
        std::vector<float> extinction;         ///< extinction coefficient (1/m)
        std::vector<float> sun_transmittance;  ///< transmittance of the dust toward the Sun, terrain ignored
        std::vector<float> sun_visibility;     ///< 1 where the terrain leaves the Sun in view, else 0
        ChVector3d sun_dir;                    ///< unit direction toward the Sun the grid was traced for
        std::uint64_t version = 0;             ///< incremented at each UpdateGrid
        size_t Index(int i, int j, int k) const { return size_t(i) + size_t(nx) * (size_t(j) + size_t(ny) * size_t(k)); }
    };

    /// Mass budget (kg) since construction.
    struct Stats {
        double emitted = 0;    ///< released
        double aloft = 0;      ///< in flight at the last Update
        double landed = 0;     ///< back on the ground
        double expired = 0;    ///< dropped after max_flight_time
        double skipped = 0;    ///< not released because max_particles were aloft
        size_t particles = 0;  ///< super-particles in flight at the last Update
    };

    ChDustField(const Params& params, HeightFunction ground);

    const Params& GetParams() const { return m_params; }

    void SetWheelEmission(const WheelEmission& emission) { m_wheel = emission; }
    const WheelEmission& GetWheelEmission() const { return m_wheel; }

    void SetGridSpec(const GridSpec& spec) { m_spec = spec; }
    const GridSpec& GetGridSpec() const { return m_spec; }

    /// Release the soil the given wheels throw up over the interval [time, time + dt), at release times spread
    /// over it. The wheels keep their state over the interval, so call this once per physics step.
    void EmitFromWheels(const std::vector<WheelState>& wheels, double time, double dt);

    /// Release one super-particle.
    void Emit(const ChVector3d& pos, const ChVector3d& vel, double mass, int bin, double time);

    /// Drop the grains that have landed or flown for too long by the given time.
    void Update(double time);

    /// Bin the grains aloft at `time` on a grid around `center` and trace it toward the Sun (unit direction).
    /// The grid snaps to whole voxels, so it does not shimmer as the center moves.
    void UpdateGrid(double time, const ChVector3d& center, const ChVector3d& sun_dir);

    const Grid& GetGrid() const { return m_grid; }
    const Stats& GetStats() const { return m_stats; }
    const std::vector<Particle>& GetParticles() const { return m_particles; }

    /// Position of a particle at a time.
    ChVector3d Position(const Particle& p, double time) const;

    /// Extinction per unit mass of a size bin (m^2/kg).
    double MassExtinction(int bin) const;

    /// Fill the grid's Sun transmittance from its extinction, for a unit direction toward the Sun. Rays from the
    /// Sun one voxel apart are carried slice by slice across the dominant axis of the direction, each summing its
    /// own optical depth, and each voxel interpolates between the four rays around it: a few lookups per voxel, and
    /// no blurring that builds up along the rays. Light enters the grid unattenuated.
    static void TraceTransmittance(Grid& grid, const ChVector3d& sun_dir);

    /// Fill the grid's Sun visibility past the terrain. Each column of voxels gets the height below which the
    /// terrain hides the Sun, S = max over d of (ground d toward the Sun - d tan(Sun elevation)), found by a
    /// sweep over the columns toward the Sun; ground beyond the grid is traced up to `shadow_distance`. Voxels
    /// fade from hidden to lit over one voxel around that height. Visibility depends only on where the grid is
    /// and on the Sun, so UpdateGrid redoes it only when either moves.
    static void TraceVisibility(Grid& grid, const ChVector3d& sun_dir, const HeightFunction& ground, double shadow_distance);

    /// Write the grid's optical depth summed along z as an 8-bit PGM image, as 1 - exp(-tau): how much of the
    /// light passing straight down through the dust it blocks.
    bool WriteOpticalDepthImage(const std::string& filename) const;

  private:
    // Ground height from a cache on a lattice of the grid's voxel size.
    double Ground(double x, double y);
    double LatticeHeight(int i, int j);

    Params m_params;
    WheelEmission m_wheel;
    GridSpec m_spec;
    HeightFunction m_ground;

    std::vector<Particle> m_particles;
    std::vector<double> m_mass_extinction;
    std::vector<double> m_bin_share;  // share of the super-particles each bin gets: its share of the cross section
    std::vector<double> m_carry;      // fractional super-particles per wheel and bin, carried to the next call
    std::mt19937_64 m_rng;
    Grid m_grid;
    Stats m_stats;

    // Ground heights on a lattice of the grid's voxel size, cached in a direct-mapped table: a lattice point has one
    // slot, which holds the last point that hashed to it
    struct HeightSlot {
        std::int32_t i = INT32_MIN;
        std::int32_t j = INT32_MIN;
        float height = 0;
    };
    double m_lattice = 0;  // spacing of the height cache (m)
    std::vector<HeightSlot> m_heights;

    // Where the grid's visibility was last traced, and toward which Sun
    ChVector3d m_vis_origin;
    ChVector3d m_vis_sun;
    bool m_vis_valid = false;
};

/// @} planet_module

}  // namespace planet
}  // namespace chrono

#endif
