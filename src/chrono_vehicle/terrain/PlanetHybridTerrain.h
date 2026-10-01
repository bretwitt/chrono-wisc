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
// Planet terrain that runs SCM soil where SCM holds and switches to CRM soil
// in a window around the rover where it does not: on steep ground and at high
// wheel slip.
//
// =============================================================================

#ifndef PLANET_HYBRID_TERRAIN_H
#define PLANET_HYBRID_TERRAIN_H

#include <memory>
#include <vector>

#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChSystem.h"

#include "chrono_vehicle/ChApiVehicle.h"
#include "chrono_vehicle/terrain/PlanetCRMWindow.h"
#include "chrono_vehicle/terrain/PlanetSCMTerrain.h"
#include "chrono_vehicle/terrain/SCMSindyResidual.h"

#include "chrono_planet/ChPlanetSurface.h"
#include "chrono_planet/ChSiteFrame.h"
#include "chrono_planet/filters/ChDeformationFilter.h"

namespace chrono {
namespace vehicle {

/// @addtogroup vehicle_terrain
/// @{

/// SCM soil (PlanetSCMTerrain) where it holds, and CRM soil (PlanetCRMWindow) in a window following the rover where it
/// does not. SCM has no soil flow: it cannot capture a wheel digging in or pushing soil aside, so it is trusted on
/// gentle ground at low slip. When the ground ahead of the rover rises steeply, or its wheels slip, for longer than a
/// dwell time, a CRM window is seeded around the rover from the ground as SCM left it (the two models share a
/// deformation filter: SCM publishes its ruts into it, and the window seeds from it). The window first runs without
/// loading the wheels (warm-up), while SCM still carries them; then the load moves from SCM to CRM over a blend time,
/// each model's forces on the wheels scaled by its share. On gentle ground at low slip the window's ruts and berms go
/// back into the filter, SCM takes them as its ground (SCMTerrain::SetModifiedNodes: the soil there starts without
/// compaction history), the load moves back over the blend time, and the window is dropped.
///
/// The SCM terrain must be initialized with the rover's wheels and set with the deformation filter, whose spacing must
/// be SCM's grid spacing. Advance the coupled system with Advance instead of stepping the multibody system.
class CH_VEHICLE_API PlanetHybridTerrain {
  public:
    /// Which soil model carries the wheels.
    enum class Mode {
        SCM,         ///< SCM alone
        CRM_WARMUP,  ///< CRM runs without load, SCM carries the wheels
        TO_CRM,      ///< the load moves from SCM to CRM
        CRM,         ///< CRM alone
        TO_SCM       ///< the load moves from CRM to SCM
    };

    /// When to switch.
    struct Policy {
        bool enabled = true;       ///< false: SCM alone throughout
        double lookahead = 1.5;    ///< distance ahead of the rover the slope is measured over (m)
        double enter_slope = 10;   ///< uphill slope ahead that calls for CRM (deg)
        double exit_slope = 6;     ///< slope ahead below which SCM may take over again (deg)
        double enter_slip = 0.4;   ///< wheel slip that calls for CRM
        double exit_slip = 0.2;    ///< wheel slip below which SCM may take over again
        double min_rim_speed = 0.02;  ///< wheel rim speed below which slip is not measured (m/s)
        double slip_filter = 0.3;  ///< time constant of the slip measurement's low-pass filter (s)
        double enter_dwell = 0.2;  ///< time a call for CRM must hold before the switch (s)
        double exit_dwell = 2.0;   ///< time the ground must stay gentle before switching back (s)
        double warmup = 0.3;       ///< time CRM runs without load before it takes any (s)
        double enter_uncertainty = 0;  ///< uncertainty of the trust model's correction that calls for CRM (0: not used)
        double exit_uncertainty = 0;   ///< uncertainty below which SCM may take over again
        double blend = 0.3;        ///< time over which the load moves between the models (s)
    };

    /// Hybrid terrain over `surface` in the `site` frame. `scm` must be initialized, with `ruts` as its deformation
    /// filter; the CRM windows have particle spacing `spacing` (m) and are configured by `setup` (see PlanetCRMWindow).
    PlanetHybridTerrain(ChSystem& sys,
                        std::shared_ptr<PlanetSCMTerrain> scm,
                        std::shared_ptr<const planet::ChPlanetSurface> surface,
                        const planet::ChSiteFrame& site,
                        std::shared_ptr<planet::ChDeformationFilter> ruts,
                        double spacing,
                        PlanetCRMWindow::Setup setup);
    ~PlanetHybridTerrain();

    void SetPolicy(const Policy& policy) { m_policy = policy; }
    const Policy& GetPolicy() const { return m_policy; }

    /// Size of the CRM window (see PlanetCRMWindow::SetWindow).
    void SetWindow(double length, double width, double depth, double margin);

    /// Raise berms beside the ruts CRM makes (see PlanetCRMWindow::SetBerms; default: false).
    void SetBerms(bool val) { m_berms = val; }

    /// Largest step of the multibody system while SCM alone carries the wheels (s, default: 1e-3).
    void SetStepSizeSCM(double step) { m_step_scm = step; }

    /// A wheel of the rover: its slip is measured, and the CRM window stamps its ruts. Wheels spin about their
    /// reference frame's y axis.
    void AddWheel(std::shared_ptr<ChBody> wheel, double radius, double width);

    /// Collect training samples for an SCM residual model: while CRM carries any of the load (after the warm-up), each
    /// Advance adds a sample per wheel near the ground, of CRM's wrench on it less SCM's, both before the
    /// force shares, with the wheel's state as SCM sees it. SCM keeps computing its forces under CRM, so the two are
    /// taken with the wheel in the same state. The model needs its wheels added and SCM's terrain set; it may also be
    /// registered with SCM (SCMTerrain::RegisterContactForceCorrection) as it collects, which does not change the
    /// samples. Null stops collecting.
    void SetResidualRecorder(std::shared_ptr<SCMSindyResidual> model) { m_recorder = std::move(model); }

    /// The residual model correcting SCM, whose trust also calls for CRM (Policy::enter_uncertainty): where its
    /// ensemble disagrees, or a wheel's state lies outside its training range, CRM is brought in (and its recordings,
    /// with SetResidualRecorder, extend the training data where it is needed). The uncertainty of a wheel is the larger
    /// of its ensemble spread and its distance outside the training range (SCMSindyResidual::GetWheelTrust).
    void SetTrustModel(std::shared_ptr<SCMSindyResidual> model) { m_trust = std::move(model); }

    /// Largest uncertainty over the wheels, as last measured (0 without a trust model).
    double GetUncertainty() const { return m_uncertainty; }

    /// Start on SCM, following `chassis`, whose x axis points forward.
    void Initialize(std::shared_ptr<ChBody> chassis);

    /// Advance the rover and the soil by `step` (s), switching models as the policy says.
    void Advance(double step);

    /// Write the ruts into the deformation filter, from the model that carries the wheels. Returns the nodes changed.
    size_t Publish();

    /// Height of the ground with the ruts (site x/y to z, m).
    double GetHeight(double x, double y) const;

    Mode GetMode() const { return m_mode; }
    static const char* GetModeName(Mode mode);

    /// Share of the wheels' load CRM carries (0 to 1).
    double GetCRMShare() const { return m_share; }

    /// Filtered wheel slip, and uphill slope of the ground ahead (deg), as last measured.
    double GetSlip() const { return m_slip; }
    double GetSlope() const { return m_slope; }

    /// The CRM window, while there is one (null on SCM alone).
    PlanetCRMWindow* GetCRMWindow() const { return m_window.get(); }

    /// Number of times CRM was brought in.
    int GetNumSwitches() const { return m_switches; }

    std::shared_ptr<PlanetSCMTerrain> GetSCM() const { return m_scm; }

  private:
    struct Wheel {
        std::shared_ptr<ChBody> body;
        double radius;
        double width;
    };

    void Measure(double step);
    void Record();
    void StartCRM();
    void GiveBackToSCM();
    void SetShare(double share);

    ChSystem& m_sys;
    std::shared_ptr<PlanetSCMTerrain> m_scm;
    std::shared_ptr<const planet::ChPlanetSurface> m_surface;
    planet::ChSiteFrame m_site;
    std::shared_ptr<planet::ChDeformationFilter> m_ruts;
    double m_spacing;
    PlanetCRMWindow::Setup m_setup;

    Policy m_policy;
    double m_length = 5.0, m_width = 3.5, m_depth = 0.25, m_margin = 1.3;
    bool m_berms = false;
    double m_step_scm = 1e-3;

    std::shared_ptr<ChBody> m_chassis;
    std::vector<Wheel> m_wheels;
    std::unique_ptr<PlanetCRMWindow> m_window;
    std::shared_ptr<SCMSindyResidual> m_recorder;
    std::shared_ptr<SCMSindyResidual> m_trust;
    double m_uncertainty = 0;

    Mode m_mode = Mode::SCM;
    double m_mode_time = 0;   // time in the current mode (s)
    double m_hold_time = 0;   // time the switch condition has held (s)
    double m_share = 0;       // CRM's share of the load
    double m_slip = 0;
    double m_slope = 0;
    int m_switches = 0;
};

/// @} vehicle_terrain

}  // namespace vehicle
}  // namespace chrono

#endif
