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
// Excavation at a lunar work site with CRM soil. As in demo_PLANET_Excavation,
// the ground of a work site near the Apollo 17 landing site is a signed
// distance field (ChSiteVolume) and a bucket digs a trench beside a berm; here
// the soil around the trench and the berm is CRM (Chrono::FSI SPH continuum
// soil, vehicle::PlanetCRMTerrain) taken from that field. The bucket moves
// through the particles, which it scoops, carries and tips out, and the soil
// resists it as the SPH solution has it. A few times per simulated second the
// particles are written back into the field, which is drawn and collided with
// as it changes. Far slower than real time: the field is where the result goes,
// for other models and for rendering.
//
// Flags: --spacing <m> (default 0.05) particle spacing, --stiffness <Pa> (default
// 1.2e5) Young's modulus of the soil, --strokes <n>
// (default 2), --time <s> ends the run, --no-vis runs without a window.
//
// =============================================================================

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iostream>

#include "chrono/assets/ChVisualShapeBox.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/utils/ChBodyGeometry.h"

#include "chrono_fsi/sph/ChFsiFluidSystemSPH.h"

#include "chrono_planet/ChPlanetSurface.h"
#include "chrono_planet/ChSiteFrame.h"
#include "chrono_planet/filters/ChDeformationFilter.h"
#include "chrono_planet/lod/ChPlanetQuadtree.h"
#include "chrono_planet/planets/moon/ChMoon.h"
#include "chrono_planet/visualization/ChPlanetVisualizationVSG.h"
#include "chrono_planet/visualization/ChSiteVolumeVisualizationVSG.h"
#include "chrono_planet/volume/ChSiteVolume.h"
#include "chrono_planet/volume/ChSiteVolumeShapes.h"

#include "chrono_vehicle/terrain/PlanetCRMTerrain.h"

#include "chrono_vsg/ChVisualSystemVSG.h"

#include "PlanetDemoSetup.h"
#include "PlanetDigCycle.h"

using namespace chrono;
using namespace chrono::planet;
using namespace chrono::vehicle;
using namespace chrono::fsi::sph;

// Landing site (Apollo 17 region), quadtree zoom of the physics surface
const double site_lon = 30.75;
const double site_lat = 20.19;
const int zoom = 15;

const double gravity = moon::kGravity;

// Work site, and the window of it that is CRM soil: the trench and the berm
const ChSiteRegion work_site(-3.5, -1.5, 1.5, 3.0);
const ChSiteRegion crm_window(-2.6, -0.7, 0.2, 2.5);
const double crm_depth = 0.35;  // soil below the lowest ground in the window (m)

// Bucket: 0.6 m long, 0.8 m wide, 0.4 m tall, open ahead and on top. Its walls are this many particle spacings
// thick, as CRM needs of the markers bounding it.
const ChVector3d bucket_size(0.6, 0.8, 0.4);
const double bucket_wall_spacings = 2.4;

// Time steps: SPH (variable, at most this), and the exchange with the multibody system
const double cfd_step = 5e-4;
const double exchange_step = 5e-3;
const double publish_period = 0.2;  // particles written into the site volume this often (s)

int main(int argc, char* argv[]) {
    std::cout << "Copyright (c) 2026 projectchrono.org\nChrono version: " << CHRONO_VERSION << std::endl;

    double spacing = 0.05;
    // Young's modulus (Pa). CRM's sound speed follows from it, and its time step from the sound speed; 1.2e5 gives
    // about 8 m/s, ten times what the bucket and the soil it pours reach, which is all SPH needs.
    double stiffness = 1.2e5;
    int strokes = 2;
    double end_time = 1e300;
    bool use_vis = true;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--spacing") && i + 1 < argc)
            spacing = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--stiffness") && i + 1 < argc)
            stiffness = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--strokes") && i + 1 < argc)
            strokes = std::max(1, std::atoi(argv[++i]));
        else if (!std::strcmp(argv[i], "--time") && i + 1 < argc)
            end_time = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--no-vis"))
            use_vis = false;
    }

    // Surface, site frame, and the work site's ground as a volume
    auto surface = CreateSiteSurface({}, zoom);
    if (!surface)
        return 1;
    ChSiteFrame site = surface->MakeSiteFrame(site_lon, site_lat);
    auto height_at = [&](double x, double y) {
        double lon, lat;
        site.ToLonLat(x, y, lon, lat);
        return surface->GetElevation(lon, lat) - site.GetOriginElevation();
    };
    ChSiteVolume::Params vparams;
    vparams.voxel = 0.05;
    vparams.depth = 1.0;
    vparams.height = 1.5;
    vparams.margin = 0.15;
    auto volume = chrono_types::make_shared<ChSiteVolume>(work_site, vparams, height_at);
    const double initial_soil = volume->GetSoilVolume();
    auto dug = chrono_types::make_shared<ChDeformationFilter>(site, vparams.voxel);
    auto world = chrono_types::make_shared<ChPlanetQuadtree>(surface->CreateView(dug), 10);

    // Multibody system: the bucket, moved along the dig cycle
    ChSystemNSC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, -gravity));
    DigPlan plan;
    plan.bucket_size = bucket_size;
    plan.trench_x0 = -2.0;
    plan.stroke_length = 0.8;
    plan.strokes_per_pass = strokes;
    plan.passes = 1;
    plan.pass_depth = 0.12;
    plan.berm_y = 1.7;
    plan.plunge_pitch = 0.5;  // edge first: a flat floor pressed into the soil stalls the solver's time step
    DigCycle cycle(volume, plan);

    auto bucket = chrono_types::make_shared<ChBody>();
    bucket->SetFixed(true);
    bucket->EnableCollision(false);
    bucket->SetPos(cycle.Start().GetPos());
    sys.AddBody(bucket);
    auto geometry = chrono_types::make_shared<utils::ChBodyGeometry>();
    geometry->materials.push_back(ChContactMaterialData());
    {
        const ChVector3d s = bucket_size;
        const double t = bucket_wall_spacings * spacing;
        geometry->coll_boxes.push_back(utils::ChBodyGeometry::BoxShape(ChVector3d(0, 0, -0.5 * (s.z() - t)), QUNIT, ChVector3d(s.x(), s.y(), t)));
        geometry->coll_boxes.push_back(utils::ChBodyGeometry::BoxShape(ChVector3d(-0.5 * (s.x() - t), 0, 0), QUNIT, ChVector3d(t, s.y(), s.z())));
        geometry->coll_boxes.push_back(utils::ChBodyGeometry::BoxShape(ChVector3d(0, -0.5 * (s.y() - t), 0), QUNIT, ChVector3d(s.x(), t, s.z())));
        geometry->coll_boxes.push_back(utils::ChBodyGeometry::BoxShape(ChVector3d(0, 0.5 * (s.y() - t), 0), QUNIT, ChVector3d(s.x(), t, s.z())));
        auto steel = chrono_types::make_shared<ChVisualMaterial>();
        steel->SetDiffuseColor(ChColor(0.75f, 0.55f, 0.15f));
        for (const auto& box : geometry->coll_boxes) {
            auto shape = chrono_types::make_shared<ChVisualShapeBox>(box.dims);
            shape->AddMaterial(steel);
            bucket->AddVisualShape(shape, ChFrame<>(box.pos, box.rot));
        }
    }

    // CRM soil over the window, from the volume
    PlanetCRMTerrain crm(sys, spacing, volume);
    crm.SetVerbose(true);
    crm.SetGravitationalAcceleration(ChVector3d(0, 0, -gravity));
    crm.SetStepSizeCFD(cfd_step);
    ChFsiFluidSystemSPH::SoilProperties soil;
    soil.density = 1600;
    soil.Young_modulus = stiffness;
    soil.Poisson_ratio = 0.3;
    soil.mu_I0 = 0.04;
    soil.mu_fric_s = 0.7;
    soil.mu_fric_2 = 0.7;
    soil.average_diam = 0.005;
    soil.cohesion_coeff = 1e3;
    crm.SetCrmSPH(soil);
    ChFsiFluidSystemSPH::SPHParameters sph;
    sph.integration_scheme = IntegrationScheme::RK2;
    sph.initial_spacing = spacing;
    sph.d0_multiplier = 1.3;
    sph.free_surface_threshold = 2.4;
    sph.artificial_viscosity = 0.5;
    sph.viscosity_method = ViscosityMethod::ARTIFICIAL_BILATERAL;
    sph.boundary_method = BoundaryMethod::ADAMI;
    sph.use_variable_time_step = true;
    sph.num_proximity_search_steps = 10;  // neighbor lists rebuilt every 10 steps
    crm.SetSPHParameters(sph);
    crm.AddRigidBody(bucket, geometry, false);
    crm.SetActiveDomain(ChVector3d(1.0, 1.2, 0.9));
    // Only a floor under the soil: away from the bucket the soil is frozen (active domains) and holds itself up, so
    // walls around the window would only add markers to every step
    crm.ConstructFromVolume(crm_window, crm_depth, BoxSide::Z_NEG);
    crm.Initialize();
    std::cout << "CRM soil: " << crm.GetNumSPHParticles() << " particles, " << crm.GetNumBoundaryBCEMarkers() << " boundary markers" << std::endl;

    // The site's ground, drawn and collided with as it changes
    auto shapes = chrono_types::make_shared<ChSiteVolumeShapes>(&sys, volume);
    shapes->SetUseVisualModel(!use_vis);
    auto regolith = chrono_types::make_shared<ChVisualMaterial>();
    regolith->SetDiffuseColor(ChColor(0.42f, 0.41f, 0.39f));
    regolith->SetRoughness(0.95f);
    auto disturbed = chrono_types::make_shared<ChVisualMaterial>();
    disturbed->SetDiffuseColor(ChColor(0.30f, 0.29f, 0.28f));
    disturbed->SetRoughness(0.95f);
    shapes->SetMaterial(regolith);
    shapes->SetDisturbedMaterial(disturbed);
    shapes->Update();

    std::shared_ptr<vsg3d::ChVisualSystemVSG> vis;
    if (use_vis) {
        auto vis_planet = chrono_types::make_shared<ChPlanetVisualizationVSG>(world, site);
        vis_planet->AddHole(volume->GetHole());
        vis = chrono_types::make_shared<vsg3d::ChVisualSystemVSG>();
        vis->AttachSystem(&sys);
        vis->AttachPlugin(vis_planet);
        vis->AttachPlugin(chrono_types::make_shared<ChSiteVolumeVisualizationVSG>(shapes));
        vis->SetWindowTitle("Excavation in CRM soil at a lunar work site");
        vis->SetWindowSize(1280, 800);
        vis->SetBackgroundColor(ChColor(0.01f, 0.01f, 0.02f));
        vis->AddCamera(ChVector3d(-3.0, -3.5, 2.5), ChVector3d(-1.2, 0.8, 0));
        vis->SetCameraVertical(CameraVerticalDir::Z);
        vis->SetLightIntensity(1.0f);
        vis->SetLightDirection(1.2, 0.5);
        vis->Initialize();
    }

    // Simulation loop: the bucket's pose for the next exchange, then the coupled step
    ChVector3d pos = cycle.Start().GetPos();
    double pitch = 0;
    double next_publish = 0, next_report = 0;
    const auto wall_start = std::chrono::steady_clock::now();
    while (sys.GetChTime() < end_time) {
        const double time = sys.GetChTime();
        if (vis && !vis->Run())
            break;

        ChVector3d velocity;
        cycle.Advance(exchange_step, pos, pitch, velocity);
        bucket->SetPos(pos);
        bucket->SetRot(QuatFromAngleY(pitch));
        bucket->SetPosDt(velocity);

        if (time >= next_publish) {
            const double change = crm.PublishToVolume();
            shapes->Update();
            volume->PublishTopSurface(*dug, site);
            if (vis)
                vis->Render();
            next_publish += publish_period;
            if (time >= next_report) {
                const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
                std::cout << "t = " << time << " s (" << time / std::max(wall, 1e-9) << "x real time), stroke "
                          << std::min(cycle.GetStroke() + 1, cycle.GetNumStrokes()) << "/" << cycle.GetNumStrokes()
                          << ", soil force on the bucket " << crm.GetFsiBodyForce(bucket).Length() << " N, soil in the site "
                          << (volume->GetSoilVolume() - initial_soil) * 1e3 << " L from the start (last write " << change * 1e3
                          << " L), trench floor " << volume->GetTopHeight(-1.6, 0) - volume->GetInitialHeight(-1.6, 0) << " m"
                          << std::endl;
                next_report += 1.0;
            }
        }

        crm.DoStepDynamics(exchange_step);
        if (cycle.Done() && !vis)
            break;
    }
    return 0;
}
