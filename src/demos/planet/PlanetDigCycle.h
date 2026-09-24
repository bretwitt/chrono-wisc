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
// A scripted dig cycle for a bucket at a planet work site: a trench cut in
// passes of strokes, each load tipped out over a berm beside it.
//
// =============================================================================

#ifndef PLANET_DIG_CYCLE_H
#define PLANET_DIG_CYCLE_H

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "chrono/core/ChFrame.h"

#include "chrono_planet/volume/ChSiteVolume.h"

/// Where and how a trench is dug: along x at y = 0, from trench_x0, in passes each pass_depth deeper than the
/// last, each pass in strokes; each load is tipped out over a berm along y = berm_y.
struct DigPlan {
    chrono::ChVector3d bucket_size = chrono::ChVector3d(0.4, 0.6, 0.3);
    double trench_x0 = -2.4;
    double stroke_length = 0.8;
    int strokes_per_pass = 4;
    int passes = 3;
    double pass_depth = 0.08;   // m
    double berm_y = 1.7;
    double drag_speed = 0.3;    // m/s through the soil
    double move_speed = 0.8;    // m/s through the air
    double tip_angle = 2.0;     // rad the bucket tips over the berm
    double tip_hold = 1.5;      // s the bucket is held tipped
    double plunge_pitch = 0;    // rad the bucket tips its cutting edge down to enter the soil (0: level)
};

/// The bucket's motion for a dig plan: legs between poses, each at a speed; a leg can follow the original ground at a
/// depth below it, and end in a dump.
class DigCycle {
  public:
    DigCycle(std::shared_ptr<chrono::planet::ChSiteVolume> volume, const DigPlan& plan) : m_volume(std::move(volume)), m_plan(plan) {
        using chrono::ChVector3d;
        auto ground = [&](double x, double y) { return m_volume->GetInitialHeight(x, y); };
        const double half = 0.5 * plan.bucket_size.z();
        for (int p = 0; p < plan.passes; ++p) {
            const double depth = plan.pass_depth * (p + 1);
            for (int s = 0; s < plan.strokes_per_pass; ++s) {
                const double x0 = plan.trench_x0 + s * plan.stroke_length, x1 = x0 + plan.stroke_length;
                const double xm = 0.5 * (x0 + x1);
                // Above the start of the stroke, down into the soil, along it, and out
                m_legs.push_back({ChVector3d(x0 - 0.3, 0, ground(x0 - 0.3, 0) + 0.6), plan.plunge_pitch, plan.move_speed});
                m_legs.push_back({ChVector3d(x0, 0, ground(x0, 0) - depth + half), plan.plunge_pitch, 0.2});
                if (plan.plunge_pitch != 0)
                    m_legs.push_back({ChVector3d(x0, 0, ground(x0, 0) - depth + half), 0, 0.5});  // level, in the soil
                m_legs.push_back({ChVector3d(x1, 0, 0), 0, plan.drag_speed, depth});
                m_legs.push_back({ChVector3d(x1, 0, ground(x1, 0) + 0.6), 0, 0.3});
                // Over the berm, tipped over and held while the load pours out, and back
                const ChVector3d over(xm, plan.berm_y, ground(xm, plan.berm_y) + 1.3);
                m_legs.push_back({over, 0, plan.move_speed});
                m_legs.push_back({over, -plan.tip_angle, 2.0, -1, true, plan.tip_hold});
                m_legs.push_back({over, 0, 2.0});
            }
        }
        const double xp = plan.trench_x0 - 1.0;
        m_legs.push_back({ChVector3d(xp, -1.0, ground(xp, -1.0) + 1.0), 0, plan.move_speed});
    }

    /// Bucket pose at the start.
    chrono::ChFrame<> Start() const { return chrono::ChFrame<>(m_legs.front().to, chrono::QUNIT); }

    /// Advance by dt from the current pose (position of the bucket center and its rotation about y); returns true
    /// when a leg ending in a dump was just completed.
    bool Advance(double dt, chrono::ChVector3d& pos, double& pitch, chrono::ChVector3d& velocity) {
        velocity = chrono::VNULL;
        if (Done())
            return false;
        const Leg& leg = m_legs[m_leg];
        const chrono::ChVector3d start = pos;
        // Holding at the end of a leg
        if (m_holding) {
            m_held += dt;
            if (m_held < leg.hold)
                return false;
            m_holding = false;
            return Finish(leg);
        }
        // Legs that only turn
        if ((leg.to - pos).Length() < 1e-9 && std::abs(leg.pitch - pitch) > 1e-9) {
            const double step = leg.speed * dt;
            pitch += std::clamp(leg.pitch - pitch, -step, step);
            if (std::abs(leg.pitch - pitch) <= 1e-9)
                return Finish(leg);
            return false;
        }
        chrono::ChVector3d target = leg.to;
        if (leg.depth >= 0)
            target.z() = m_volume->GetInitialHeight(target.x(), target.y()) - leg.depth + 0.5 * m_plan.bucket_size.z();
        const chrono::ChVector3d d = target - pos;
        const double dist = d.Length();
        const double step = leg.speed * dt;
        if (dist <= step) {
            pos = target;
        } else {
            pos += d * (step / dist);
            if (leg.depth >= 0)
                pos.z() = m_volume->GetInitialHeight(pos.x(), pos.y()) - leg.depth + 0.5 * m_plan.bucket_size.z();
        }
        velocity = (pos - start) / dt;
        if ((pos - target).Length() < 1e-9 && std::abs(leg.pitch - pitch) <= 1e-9)
            return Finish(leg);
        return false;
    }

    bool Done() const { return m_leg >= m_legs.size(); }
    int GetStroke() const { return static_cast<int>(m_leg / (m_plan.plunge_pitch != 0 ? kLegsPerStroke + 1 : kLegsPerStroke)); }
    int GetNumStrokes() const { return m_plan.passes * m_plan.strokes_per_pass; }

  private:
    static constexpr size_t kLegsPerStroke = 7;
    struct Leg {
        chrono::ChVector3d to;  // bucket center at the end of the leg (z ignored when following the ground)
        double pitch;           // rotation about y at the end of the leg (rad)
        double speed;           // m/s, or rad/s for legs that only turn
        double depth = -1;      // follow the original ground with the bucket bottom this far below it, if >= 0
        bool dump = false;      // let what is left in the bucket fall at the end of the leg
        double hold = 0;        // stay at the end of the leg this long (s)
    };

    bool Finish(const Leg& leg) {
        if (leg.hold > 0 && !m_holding && m_held == 0) {
            m_holding = true;
            return false;
        }
        m_held = 0;
        ++m_leg;
        return leg.dump;
    }

    std::shared_ptr<chrono::planet::ChSiteVolume> m_volume;
    DigPlan m_plan;
    std::vector<Leg> m_legs;
    size_t m_leg = 0;
    bool m_holding = false;
    double m_held = 0;
};

#endif
