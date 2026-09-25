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
// Tests of ChDustField: ballistic flight, the mass budget of wheel emission,
// the direction wheels throw soil, and the voxel grid and its Sun tracing.
//
// =============================================================================

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

#include "gtest/gtest.h"

#include "chrono/utils/ChConstants.h"

#include "chrono_vehicle/terrain/ChDustField.h"

using namespace chrono;
using namespace chrono::vehicle;

namespace {

const double kGravity = 1.62;

double Flat(double, double) {
    return 0.0;
}

ChDustField::Params TestParams() {
    ChDustField::Params params;
    params.gravity = kGravity;
    return params;
}

// A wheel of radius r resting on flat ground, axle along +y, spinning at omega and moving at v along +x
ChDustField::WheelState Wheel(double r, double omega, double v) {
    return {ChVector3d(0, 0, r), ChVector3d(0, 1, 0), ChVector3d(v, 0, 0), omega, r, 0.3};
}

// Distance from p along the unit direction L to the boundary of the box [lo, hi]
double ExitDistance(const ChVector3d& p, const ChVector3d& L, const ChVector3d& lo, const ChVector3d& hi) {
    double t = std::numeric_limits<double>::max();
    for (int a = 0; a < 3; ++a) {
        if (std::abs(L[a]) < 1e-12)
            continue;
        const double bound = L[a] > 0 ? hi[a] : lo[a];
        t = std::min(t, (bound - p[a]) / L[a]);
    }
    return t;
}

}  // namespace

// A grain lands after 2 v_z / g, at the top of its arc halfway there.
TEST(ChDustField, BallisticFlight) {
    ChDustField dust(TestParams(), Flat);
    const ChVector3d v(2.0, 0.5, 2.0);
    dust.Emit(ChVector3d(0, 0, 0), v, 1.0, 0, 0.0);
    const double flight = 2 * v.z() / kGravity;

    const ChVector3d apex = dust.Position(dust.GetParticles()[0], 0.5 * flight);
    EXPECT_NEAR(apex.z(), v.z() * v.z() / (2 * kGravity), 1e-12);
    EXPECT_NEAR(apex.x(), v.x() * 0.5 * flight, 1e-12);

    dust.Update(flight - 0.01);
    EXPECT_EQ(dust.GetStats().particles, 1u);
    dust.Update(flight + 0.01);
    EXPECT_EQ(dust.GetStats().particles, 0u);
    EXPECT_DOUBLE_EQ(dust.GetStats().landed, 1.0);
}

// A grain not yet released is neither moved nor dropped.
TEST(ChDustField, FutureRelease) {
    ChDustField dust(TestParams(), Flat);
    dust.Emit(ChVector3d(0, 0, 0.1), ChVector3d(0, 0, -1), 1.0, 0, 5.0);
    dust.Update(1.0);
    EXPECT_EQ(dust.GetStats().particles, 1u);
}

// Every kilogram released is aloft, landed, expired or skipped, and a wheel spinning in place releases
// bulk_density * width * loose_depth * (1 + slip_gain) * rim speed per second.
TEST(ChDustField, WheelMassBudget) {
    auto params = TestParams();
    params.max_flight_time = 0.5;
    ChDustField dust(params, Flat);
    const double r = 0.25, omega = 8.0, dt = 1e-3;
    const std::vector<ChDustField::WheelState> wheels = {Wheel(r, omega, 0.0)};
    double time = 0;
    for (int n = 0; n < 1000; ++n) {
        dust.EmitFromWheels(wheels, time, dt);
        time += dt;
        dust.Update(time);
    }
    const auto& s = dust.GetStats();
    EXPECT_GT(s.landed, 0.0);
    EXPECT_NEAR(s.emitted, s.aloft + s.landed + s.expired + s.skipped, 1e-9 * s.emitted);

    const auto& e = dust.GetWheelEmission();
    const double rate = e.bulk_density * 0.3 * e.loose_depth * (1 + e.slip_gain) * omega * r;
    EXPECT_NEAR(s.emitted, rate * time, 0.01 * rate * time);
}

// A wheel spinning in place throws soil backward and up from behind it: the rim rises out of the ground there.
TEST(ChDustField, WheelThrowsBackward) {
    ChDustField dust(TestParams(), Flat);
    auto emission = dust.GetWheelEmission();
    emission.spread = 0;
    dust.SetWheelEmission(emission);
    const double r = 0.25, omega = 8.0;
    dust.EmitFromWheels({Wheel(r, omega, 0.0)}, 0.0, 0.01);
    ASSERT_FALSE(dust.GetParticles().empty());
    for (const auto& p : dust.GetParticles()) {
        EXPECT_LE(p.pos0.x(), 1e-12);
        EXPECT_GE(p.pos0.z(), -1e-12);
        EXPECT_LE(p.vel0.x(), 1e-12);
        EXPECT_GE(p.vel0.z(), -1e-12);
        EXPECT_LE(p.vel0.Length(), omega * r * (1 + 1e-12));
    }
}

// A wheel off the ground or not turning throws nothing.
TEST(ChDustField, WheelOffGround) {
    ChDustField dust(TestParams(), Flat);
    auto lifted = Wheel(0.25, 8.0, 0.0);
    lifted.center.z() += 0.1;
    dust.EmitFromWheels({lifted, Wheel(0.25, 0.0, 0.0)}, 0.0, 0.1);
    EXPECT_TRUE(dust.GetParticles().empty());
}

// Binning keeps the extinction: its integral over the grid is the particles' total cross section.
TEST(ChDustField, GridConservesExtinction) {
    auto params = TestParams();
    params.gravity = 0;
    ChDustField dust(params, [](double, double) { return -10.0; });
    ChDustField::GridSpec spec;
    spec.voxel = 0.1;
    spec.nx = spec.ny = spec.nz = 20;
    spec.below = 1.0;
    dust.SetGridSpec(spec);

    std::mt19937 rng(7);
    std::uniform_real_distribution<double> pos(-0.8, 0.8);
    double cross_section = 0;
    for (int n = 0; n < 500; ++n) {
        const int bin = n % 3;
        dust.Emit(ChVector3d(pos(rng), pos(rng), pos(rng)), ChVector3d(0, 0, 0), 1e-6 * (1 + n % 5), bin, 0.0);
        cross_section += 1e-6 * (1 + n % 5) * dust.MassExtinction(bin);
    }
    dust.UpdateGrid(0.0, ChVector3d(0, 0, 0), ChVector3d(0.3, 0.2, 0.9));
    const auto& g = dust.GetGrid();
    double sum = 0;
    for (float s : g.extinction)
        sum += s;
    EXPECT_NEAR(sum * g.voxel * g.voxel * g.voxel, cross_section, 1e-5 * cross_section);
}

// Behind a smooth cloud of dust the traced optical depth matches the extinction integrated along the ray to the
// Sun, to the accuracy of sampling the cloud at the voxel centers.
TEST(ChDustField, TransmittanceCloud) {
    ChDustField::Grid g;
    g.origin = ChVector3d(-2, -2, 0);
    g.voxel = 0.1;
    g.nx = g.ny = 40;
    g.nz = 30;
    const size_t count = size_t(g.nx) * g.ny * g.nz;
    const ChVector3d center(0.2, -0.3, 1.2);
    const double width = 0.4, peak = 3.0;
    auto sigma = [&](const ChVector3d& p) { return peak * std::exp(-(p - center).Length2() / (2 * width * width)); };
    g.extinction.resize(count);
    for (int k = 0; k < g.nz; ++k)
        for (int j = 0; j < g.ny; ++j)
            for (int i = 0; i < g.nx; ++i)
                g.extinction[g.Index(i, j, k)] = float(sigma(g.origin + ChVector3d(i + 0.5, j + 0.5, k + 0.5) * g.voxel));
    g.sun_transmittance.resize(count);
    const ChVector3d L = ChVector3d(1.0, -0.3, 0.5).GetNormalized();
    ChDustField::TraceTransmittance(g, L);

    const ChVector3d hi = g.origin + ChVector3d(g.nx, g.ny, g.nz) * g.voxel;
    double max_depth = 0;
    for (int k = 0; k < g.nz; ++k)
        for (int j = 0; j < g.ny; ++j)
            for (int i = 0; i < g.nx; ++i) {
                const ChVector3d p = g.origin + ChVector3d(i + 0.5, j + 0.5, k + 0.5) * g.voxel;
                const double d = ExitDistance(p, L, g.origin, hi);
                double exact = 0;
                const int steps = 400;
                for (int n = 0; n < steps; ++n)
                    exact += sigma(p + L * ((n + 0.5) * d / steps)) * d / steps;
                const double depth = -std::log(g.sun_transmittance[g.Index(i, j, k)]);
                max_depth = std::max(max_depth, exact);
                ASSERT_NEAR(depth, exact, 0.03 * exact + 0.02) << i << " " << j << " " << k;
            }
    // The cloud is thick enough to shade what is behind it
    EXPECT_GT(max_depth, 2.0);
}

// In dust filling the grid the traced optical depth is sigma times the distance to where the light enters the grid,
// to within a sweep step: a ray entering through a side starts at the first slice it reaches.
TEST(ChDustField, TransmittanceUniform) {
    ChDustField::Grid g;
    g.origin = ChVector3d(-1, -1, 0);
    g.voxel = 0.1;
    g.nx = g.ny = g.nz = 20;
    const size_t count = size_t(g.nx) * g.ny * g.nz;
    const double sigma = 2.0;
    g.extinction.assign(count, float(sigma));
    g.sun_transmittance.resize(count);
    const ChVector3d L = ChVector3d(1.0, -0.3, 0.5).GetNormalized();
    ChDustField::TraceTransmittance(g, L);

    const double step = g.voxel / std::abs(L.x());
    const ChVector3d hi = g.origin + ChVector3d(g.nx, g.ny, g.nz) * g.voxel;
    for (int k = 0; k < g.nz; ++k)
        for (int j = 0; j < g.ny; ++j)
            for (int i = 0; i < g.nx; ++i) {
                const ChVector3d p = g.origin + ChVector3d(i + 0.5, j + 0.5, k + 0.5) * g.voxel;
                const double d = ExitDistance(p, L, g.origin, hi);
                const double depth = -std::log(g.sun_transmittance[g.Index(i, j, k)]);
                ASSERT_NEAR(depth, sigma * d, sigma * step) << i << " " << j << " " << k;
            }
}

// A wall hides the Sun below the line grazing its top, whether the wall is inside the grid or past it.
TEST(ChDustField, VisibilityBehindWall) {
    const ChVector3d L(std::cos(CH_PI / 6), 0, std::sin(CH_PI / 6));  // Sun 30 deg up in the east
    const double tan_elev = std::tan(CH_PI / 6);
    for (double wall_x : {0.5, 1.5}) {
        ChDustField::Grid g;
        g.origin = ChVector3d(-1, -1, 0);
        g.voxel = 0.1;
        g.nx = g.ny = 20;
        g.nz = 30;
        g.sun_visibility.resize(size_t(g.nx) * g.ny * g.nz);
        ChDustField::TraceVisibility(g, L, [&](double x, double) { return x > wall_x ? 1.0 : 0.0; }, 10.0);

        // Column at x = -0.45: the first column on the wall is at x = wall_x + 0.05 inside the grid (the grid
        // samples the ground at column centers) and the wall face at wall_x past it
        const int i = 5;
        const double x = g.origin.x() + (i + 0.5) * g.voxel;
        const double face = wall_x < 1.0 ? wall_x + 0.05 : wall_x;
        const double shade = 1.0 - tan_elev * (face - x);
        for (int k = 0; k < g.nz; ++k) {
            const double z = g.origin.z() + (k + 0.5) * g.voxel;
            const float v = g.sun_visibility[g.Index(i, 10, k)];
            if (z < shade - 0.1)
                EXPECT_EQ(v, 0.f) << "wall " << wall_x << " z " << z;
            else if (z > shade + 0.1)
                EXPECT_EQ(v, 1.f) << "wall " << wall_x << " z " << z;
        }
    }
}

// The Sun below the horizon lights nothing.
TEST(ChDustField, VisibilitySunDown) {
    ChDustField::Grid g;
    g.origin = ChVector3d(0, 0, 0);
    g.voxel = 0.1;
    g.nx = g.ny = g.nz = 4;
    g.sun_visibility.assign(64, 1.f);
    ChDustField::TraceVisibility(g, ChVector3d(1, 0, -0.1).GetNormalized(), Flat, 10.0);
    for (float v : g.sun_visibility)
        EXPECT_EQ(v, 0.f);
}
