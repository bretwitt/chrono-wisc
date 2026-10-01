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

#include "chrono_vehicle/terrain/PlanetHybridTerrain.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace chrono {
namespace vehicle {

PlanetHybridTerrain::PlanetHybridTerrain(ChSystem& sys,
                                         std::shared_ptr<PlanetSCMTerrain> scm,
                                         std::shared_ptr<const planet::ChPlanetSurface> surface,
                                         const planet::ChSiteFrame& site,
                                         std::shared_ptr<planet::ChDeformationFilter> ruts,
                                         double spacing,
                                         PlanetCRMWindow::Setup setup)
    : m_sys(sys), m_scm(std::move(scm)), m_surface(std::move(surface)), m_site(site), m_ruts(std::move(ruts)), m_spacing(spacing), m_setup(std::move(setup)) {
    if (!m_scm || !m_surface || !m_ruts || !m_setup)
        throw std::invalid_argument("PlanetHybridTerrain: null SCM terrain, surface, ruts or setup");
    if (std::abs(m_ruts->GetSpacing() - m_scm->GetParams().delta) > 1e-9)
        throw std::invalid_argument("PlanetHybridTerrain: the deformation filter's spacing must be SCM's grid spacing");
}

PlanetHybridTerrain::~PlanetHybridTerrain() {}

const char* PlanetHybridTerrain::GetModeName(Mode mode) {
    switch (mode) {
        case Mode::SCM:
            return "SCM";
        case Mode::CRM_WARMUP:
            return "CRM warm-up";
        case Mode::TO_CRM:
            return "SCM to CRM";
        case Mode::CRM:
            return "CRM";
        case Mode::TO_SCM:
            return "CRM to SCM";
    }
    return "";
}

void PlanetHybridTerrain::SetWindow(double length, double width, double depth, double margin) {
    m_length = length;
    m_width = width;
    m_depth = depth;
    m_margin = margin;
}

void PlanetHybridTerrain::AddWheel(std::shared_ptr<ChBody> wheel, double radius, double width) {
    if (!wheel || !(radius > 0) || !(width > 0))
        throw std::invalid_argument("PlanetHybridTerrain::AddWheel: null wheel or bad size");
    m_wheels.push_back({std::move(wheel), radius, width});
}

void PlanetHybridTerrain::Initialize(std::shared_ptr<ChBody> chassis) {
    if (!chassis)
        throw std::invalid_argument("PlanetHybridTerrain::Initialize: null chassis");
    if (m_wheels.empty())
        throw std::invalid_argument("PlanetHybridTerrain::Initialize: no wheels");
    m_chassis = std::move(chassis);
    SetShare(0);
}

double PlanetHybridTerrain::GetHeight(double x, double y) const {
    double lon, lat;
    m_site.ToLonLat(x, y, lon, lat);
    return m_surface->GetElevation(lon, lat) - m_site.GetOriginElevation() + m_ruts->GetDelta(lon, lat, 0.0);
}

void PlanetHybridTerrain::Measure(double step) {
    // Slip: how much of the wheels' rim speed does not carry the rover forward, low-pass filtered
    ChVector3d fwd = m_chassis->GetFrameRefToAbs().GetRotMat().GetAxisX();
    double rim = 0;
    for (const auto& w : m_wheels) {
        const ChVector3d axle = w.body->GetFrameRefToAbs().GetRotMat().GetAxisY();
        rim += std::abs(w.body->GetAngVelParent().Dot(axle)) * w.radius / m_wheels.size();
    }
    const double advance = m_chassis->GetPosDt().Dot(fwd);
    const double slip = rim > m_policy.min_rim_speed ? std::clamp(1 - advance / rim, 0.0, 1.0) : 0.0;
    m_slip += (1 - std::exp(-step / m_policy.slip_filter)) * (slip - m_slip);

    // Trust in the residual correction: the least trusted wheel
    m_uncertainty = 0;
    if (m_trust) {
        for (const auto& w : m_wheels) {
            double spread, outside;
            if (m_trust->GetWheelTrust(*w.body, spread, outside))
                m_uncertainty = std::max({m_uncertainty, spread, outside});
        }
    }

    // Slope: rise of the ground from under the rover to a lookahead distance ahead of it, along its heading
    fwd.z() = 0;
    if (fwd.Length() < 1e-6)
        return;
    fwd.Normalize();
    const ChVector3d p = m_chassis->GetPos();
    const ChVector3d q = p + fwd * m_policy.lookahead;
    m_slope = std::atan2(GetHeight(q.x(), q.y()) - GetHeight(p.x(), p.y()), m_policy.lookahead) * CH_RAD_TO_DEG;
}

void PlanetHybridTerrain::SetShare(double share) {
    m_share = std::clamp(share, 0.0, 1.0);
    m_scm->SetForceScale(1 - m_share);
    if (m_window)
        m_window->SetForceScale(m_share);
}

void PlanetHybridTerrain::StartCRM() {
    // The window seeds from the ground as SCM left it
    m_scm->PublishDeformation();
    m_window = std::make_unique<PlanetCRMWindow>(m_sys, m_surface, m_site, m_ruts, m_spacing, m_setup);
    m_window->SetWindow(m_length, m_width, m_depth, m_margin);
    m_window->SetBerms(m_berms);
    for (const auto& w : m_wheels)
        m_window->AddWheel(w.body, w.radius, w.width);
    m_window->SetForceScale(0);
    m_window->Initialize(m_chassis);
    ++m_switches;
}

void PlanetHybridTerrain::GiveBackToSCM() {
    // The window's ruts and berms into the filter, and SCM takes them, over the window, as its ground
    m_window->Publish();
    const planet::ChSiteRegion& region = m_window->GetWindow();
    const double h = m_ruts->GetSpacing();
    const int i0 = static_cast<int>(std::ceil(region.min_x / h)), i1 = static_cast<int>(std::floor(region.max_x / h));
    const int j0 = static_cast<int>(std::ceil(region.min_y / h)), j1 = static_cast<int>(std::floor(region.max_y / h));
    std::vector<ChVector2i> locs;
    for (int j = j0; j <= j1; ++j)
        for (int i = i0; i <= i1; ++i)
            locs.push_back(ChVector2i(i, j));
    std::vector<double> ground;
    m_scm->GetHeightFunctor()->GetInitHeights(locs, h, ground);
    std::vector<SCMTerrain::NodeLevel> nodes;
    nodes.reserve(locs.size());
    for (size_t n = 0; n < locs.size(); ++n)
        nodes.push_back({locs[n], ground[n] + m_ruts->GetDelta(locs[n])});
    m_scm->SetModifiedNodes(nodes);
}

void PlanetHybridTerrain::Advance(double step) {
    Measure(step);
    const Policy& p = m_policy;
    m_mode_time += step;
    auto enter = [&](Mode mode) {
        m_mode = mode;
        m_mode_time = 0;
        m_hold_time = 0;
    };
    switch (m_mode) {
        case Mode::SCM:
            m_hold_time = p.enabled && (m_slope > p.enter_slope || m_slip > p.enter_slip ||
                                        (m_trust && p.enter_uncertainty > 0 && m_uncertainty > p.enter_uncertainty))
                              ? m_hold_time + step
                              : 0;
            if (m_hold_time >= p.enter_dwell) {
                StartCRM();
                enter(Mode::CRM_WARMUP);
            }
            break;
        case Mode::CRM_WARMUP:
            if (m_mode_time >= p.warmup)
                enter(Mode::TO_CRM);
            break;
        case Mode::TO_CRM:
            SetShare(m_mode_time / p.blend);
            if (m_mode_time >= p.blend)
                enter(Mode::CRM);
            break;
        case Mode::CRM:
            m_hold_time = m_slope < p.exit_slope && m_slip < p.exit_slip &&
                                  (!m_trust || p.enter_uncertainty <= 0 || m_uncertainty < p.exit_uncertainty)
                              ? m_hold_time + step
                              : 0;
            if (m_hold_time >= p.exit_dwell) {
                GiveBackToSCM();
                enter(Mode::TO_SCM);
            }
            break;
        case Mode::TO_SCM:
            SetShare(1 - m_mode_time / p.blend);
            if (m_mode_time >= p.blend) {
                m_window.reset();
                SetShare(0);
                m_scm->PublishDeformation();
                enter(Mode::SCM);
            }
            break;
    }

    if (m_window) {
        m_window->Advance(step);
        if (m_recorder && m_mode != Mode::CRM_WARMUP)
            Record();
    } else {
        double t = 0;
        while (t < step - 1e-12) {
            const double h = std::min(m_step_scm, step - t);
            m_sys.DoStepDynamics(h);
            t += h;
        }
    }
}

void PlanetHybridTerrain::Record() {
    const PlanetCRMTerrain& crm = m_window->GetTerrain();
    for (const auto& w : m_wheels) {
        ChVector3d scm_force, scm_torque, crm_force, crm_torque;
        m_scm->GetContactForceBody(w.body, scm_force, scm_torque);  // zero where SCM does not touch the wheel
        if (crm.GetSoilWrench(*w.body, crm_force, crm_torque))
            m_recorder->AddSample(*w.body, scm_force, scm_torque, crm_force, crm_torque);
    }
}

size_t PlanetHybridTerrain::Publish() {
    if (m_window)
        return m_window->Publish();
    return m_scm->PublishDeformation();
}

}  // namespace vehicle
}  // namespace chrono
