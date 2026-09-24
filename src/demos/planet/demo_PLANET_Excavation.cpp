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
// Excavation at a lunar work site. The ground over an 8 x 6 m rectangle near
// the Apollo 17 landing site is taken over by a ChSiteVolume, a signed distance
// field that can hold what digging makes and a height field cannot: a trench
// with vertical walls, and a berm of spoil beside it. A bucket, moved along a
// scripted dig cycle, cuts the trench in passes of increasing depth: each
// stroke drags it through the soil, which it carries (ChExcavationTool), lifts
// it out, and tips it over the berm. The soil stays in view the whole way
// (ChSoilParticles): cut, it fills the bucket as particles drawn as one smooth
// mass; tipped past the angle of repose, it pours out, falls, and lands, and
// becomes part of the berm where it lands. The soil resists the bucket as the
// fundamental equation of earthmoving has it, and the run reports the forces,
// the volumes cut and dumped, and how well the soil in the site balances.
//
// The quadtree terrain is drawn around the site with a hole cut where the
// volume takes over (ChPlanetVisualizationVSG::AddHole), and the volume's
// ground is drawn brick by brick as it changes (ChSiteVolumeVisualizationVSG).
// The top of the dug ground also goes into a deformation filter on the drawn
// surface, for tiles far enough away to be drawn coarsely.
//
// With --rover, a VIPER rover waits beside the site until the digging is done
// and then drives over the berm: it rides the volume's collision meshes
// (ChSiteVolumeShapes) inside the site and a planet terrain patch, with the site
// cut out of it, around it.
//
// Flags: --passes <n> (default 3) trench passes, --voxel <m> (default 0.05),
// --rover, --time <s> ends the run, --no-vis runs without a window,
// --snapshots saves the window every 5 simulated seconds (--snapshot-period <s>)
// to the demo output directory, --close puts the camera by the trench.
//
// =============================================================================

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "chrono/assets/ChVisualShapeBox.h"
#include "chrono/core/ChDataPath.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/solver/ChIterativeSolver.h"
#include "chrono/utils/ChOpenMP.h"

#include "chrono_models/robot/viper/Viper.h"

#include "chrono_planet/ChPlanetSurface.h"
#include "chrono_planet/ChSiteFrame.h"
#include "chrono_planet/filters/ChDeformationFilter.h"
#include "chrono_planet/lod/ChPlanetQuadtree.h"
#include "chrono_planet/planets/moon/ChMoon.h"
#include "chrono_planet/visualization/ChPlanetVisualizationVSG.h"
#include "chrono_planet/visualization/ChSiteVolumeVisualizationVSG.h"
#include "chrono_planet/visualization/ChSoilParticlesVisualizationVSG.h"
#include "chrono_planet/volume/ChExcavationTool.h"
#include "chrono_planet/volume/ChSiteVolume.h"
#include "chrono_planet/volume/ChSiteVolumeShapes.h"
#include "chrono_planet/volume/ChSoilParticles.h"

#include "chrono_vehicle/terrain/PlanetTerrain.h"

#include "chrono_vsg/ChGuiComponentVSG.h"
#include "chrono_vsg/ChVisualSystemVSG.h"

#include "PlanetDemoSetup.h"
#include "PlanetDigCycle.h"

using namespace chrono;
using namespace chrono::planet;
using namespace chrono::vehicle;
using namespace chrono::viper;

// Landing site (Apollo 17 region) and the quadtree zoom the physics surface is sampled at
const double site_lon = 30.75;
const double site_lat = 20.19;
const int zoom = 15;
const double root_tile_deg = 16.0;
const int view_range_tiles = 10;

// Physics
const double gravity = moon::kGravity;
const double step_size = 2e-3;
const int render_steps = 17;  // physics steps per rendered frame, about 30 frames per simulated second

// Work site: the ground the volume takes over, in site coordinates (m)
const ChSiteRegion work_site(-4.0, -3.0, 4.0, 3.0);

// Bucket: a scoop 0.6 m wide, 0.4 m long and 0.3 m tall, open ahead (+x) and on top
const ChVector3d bucket_size(0.4, 0.6, 0.3);
const double bucket_wall = 0.02;
const double bucket_capacity = 0.06;  // m^3, bank volume
const double particle_spacing = 0.06;  // soil in the bucket and in the air, as particles this far apart (m)

// Dig plan: a trench along x at y = 0, cut in strokes, the spoil dumped on a berm along y = 1.7. A stroke and
// the bucket's own length at one pass depth hold about a bucketful, 0.6 x 0.08 x 1.2 m.
DigPlan MakePlan(int passes) {
    DigPlan plan;
    plan.bucket_size = bucket_size;
    plan.trench_x0 = -2.4;
    plan.stroke_length = 0.8;
    plan.strokes_per_pass = 4;
    plan.passes = passes;
    plan.pass_depth = 0.08;
    plan.berm_y = 1.7;
    return plan;
}

// Rover
const double wheel_radius = 0.25;
const double rover_speed = 1.2;  // wheel angular speed, rad/s

// A driver that keeps the rover braked until a start time, then ramps the wheels up to a speed
class DelayedSpeedDriver : public ViperDriver {
  public:
    DelayedSpeedDriver(double speed) : m_speed(speed) {}
    void Start(double time) { m_start = time; }
    bool Started() const { return m_start < 1e300; }

  private:
    virtual DriveMotorType GetDriveMotorType() const override { return DriveMotorType::SPEED; }
    virtual void Update(double time) override {
        const double ramp = std::clamp((time - m_start) / 2.0, 0.0, 1.0);
        for (int i = 0; i < 4; ++i) {
            drive_speeds[i] = m_speed * ramp;
            steer_angles[i] = 0;
            lift_angles[i] = 0;
        }
    }
    double m_speed;
    double m_start = 1e300;
};

// -----------------------------------------------------------------------------
// GUI panel
// -----------------------------------------------------------------------------
struct SiteStats {
    int stroke = 0, strokes = 0;
    double depth = 0, force = 0, power = 0, payload = 0;
    size_t carried = 0, falling = 0, resting = 0;
    double excavated = 0, dumped = 0, spilled = 0, balance = 0;
    size_t bricks = 0, triangles = 0, rebuilt = 0;
    bool rover = false;
    double rover_x = 0;
};

class SiteStatsVSG : public vsg3d::ChGuiComponentVSG {
  public:
    SiteStatsVSG(const SiteStats& stats) : m_stats(stats) {}
    virtual void render(vsg::CommandBuffer& cb) override {
        ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f));
        ImGui::Begin("Excavation");
        ImGui::Text("Stroke %d of %d", std::min(m_stats.stroke + 1, m_stats.strokes), m_stats.strokes);
        ImGui::Text("Cutting depth: %.3f m", m_stats.depth);
        ImGui::Text("Soil resistance: %.0f N, power %.0f W", m_stats.force, m_stats.power);
        ImGui::Text("Payload: %.1f L", m_stats.payload * 1e3);
        ImGui::Text("Particles: %zu carried, %zu falling, %zu landed", m_stats.carried, m_stats.falling, m_stats.resting);
        ImGui::Separator();
        ImGui::Text("Cut: %.1f L, dumped %.1f L, pushed ahead %.1f L", m_stats.excavated * 1e3, m_stats.dumped * 1e3, m_stats.spilled * 1e3);
        ImGui::Text("Soil balance: %.2f L", m_stats.balance * 1e3);
        ImGui::Separator();
        ImGui::Text("Ground: %zu bricks, %zu triangles", m_stats.bricks, m_stats.triangles);
        ImGui::Text("Bricks rebuilt last frame: %zu", m_stats.rebuilt);
        if (m_stats.rover)
            ImGui::Text("Rover at x = %.2f m", m_stats.rover_x);
        ImGui::End();
    }

  private:
    const SiteStats& m_stats;
};

std::shared_ptr<ChVisualMaterial> Material(float r, float g, float b) {
    auto mat = chrono_types::make_shared<ChVisualMaterial>();
    mat->SetDiffuseColor(ChColor(r, g, b));
    mat->SetRoughness(0.95f);
    mat->SetMetallic(0.0f);
    return mat;
}

int main(int argc, char* argv[]) {
    std::cout << "Copyright (c) 2026 projectchrono.org\nChrono version: " << CHRONO_VERSION << std::endl;

    int passes = 3;
    double voxel = 0.05;
    bool use_rover = false;
    bool use_vis = true;
    bool snapshots = false;
    double snapshot_period = 5.0;
    bool close_camera = false;
    double end_time = 1e300;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--passes") && i + 1 < argc)
            passes = std::max(1, std::atoi(argv[++i]));
        else if (!std::strcmp(argv[i], "--voxel") && i + 1 < argc)
            voxel = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--rover"))
            use_rover = true;
        else if (!std::strcmp(argv[i], "--time") && i + 1 < argc)
            end_time = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--no-vis"))
            use_vis = false;
        else if (!std::strcmp(argv[i], "--snapshots"))
            snapshots = true;
        else if (!std::strcmp(argv[i], "--snapshot-period") && i + 1 < argc)
            snapshot_period = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--close"))
            close_camera = true;
    }

    // Surface, site frame, and the height of the undisturbed ground
    auto surface = CreateSiteSurface({}, zoom, root_tile_deg);
    if (!surface)
        return 1;
    ChSiteFrame site = surface->MakeSiteFrame(site_lon, site_lat);
    auto height_at = [&](double x, double y) {
        double lon, lat;
        site.ToLonLat(x, y, lon, lat);
        return surface->GetElevation(lon, lat) - site.GetOriginElevation();
    };

    // The drawn terrain: the planet surface with the top of the dug ground on it, for tiles drawn too coarsely to
    // be cut around the site
    auto dug = chrono_types::make_shared<ChDeformationFilter>(site, voxel);
    auto drawn = surface->CreateView(dug);
    auto world = chrono_types::make_shared<ChPlanetQuadtree>(drawn, view_range_tiles);

    // The work site's ground as a volume, taken from the ground as the finest tiles draw it, so the site and the
    // tiles around it meet at the same height
    ChSiteVolume::Params vparams;
    vparams.voxel = voxel;
    vparams.depth = 1.0;
    vparams.height = 1.5;
    vparams.margin = 0.15;
    const int drawn_zoom = world->GetMaxZoom();
    auto drawn_height_at = [&](double x, double y) {
        double lon, lat;
        site.ToLonLat(x, y, lon, lat);
        return drawn->GetElevationAtZoom(lon, lat, drawn_zoom) - site.GetOriginElevation();
    };
    const auto init_start = std::chrono::steady_clock::now();
    auto volume = chrono_types::make_shared<ChSiteVolume>(work_site, vparams, drawn_height_at);
    std::cout << "Site volume: " << volume->GetNumAllocatedBricks() << " bricks of ground at " << voxel << " m, built in "
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - init_start).count() << " s" << std::endl;
    const double initial_soil = volume->GetSoilVolume();

    // Multibody system
    ChSystemNSC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, -gravity));
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.SetSolverType(ChSolver::Type::BARZILAIBORWEIN);
    sys.GetSolver()->AsIterative()->SetMaxIterations(50);
    sys.SetNumThreads(std::min(8, ChOMP::GetNumProcs()), 1, 4);

    auto contact_mat = chrono_types::make_shared<ChContactMaterialNSC>();
    contact_mat->SetFriction(0.8f);
    contact_mat->SetRestitution(0.01f);

    // The volume's ground: drawn by the plugin below, and collided with
    auto shapes = chrono_types::make_shared<ChSiteVolumeShapes>(&sys, volume);
    shapes->SetUseVisualModel(false);
    shapes->SetMaterial(Material(0.42f, 0.41f, 0.39f));
    shapes->SetDisturbedMaterial(Material(0.35f, 0.34f, 0.32f));
    shapes->EnableCollision(contact_mat);
    shapes->Update();
    std::cout << "Site ground: " << shapes->GetNumShapes() << " brick meshes, " << shapes->GetNumTriangles() << " triangles" << std::endl;

    // Rigid terrain around the site for the rover, with the site cut out
    std::unique_ptr<PlanetTerrain> terrain;
    std::unique_ptr<Viper> viper;
    std::shared_ptr<DelayedSpeedDriver> driver;
    if (use_rover) {
        terrain = std::make_unique<PlanetTerrain>(&sys, surface, site);
        terrain->SetPatchSize(40.0);
        terrain->SetContactMaterial(contact_mat);
        terrain->AddHole(volume->GetHole(), [volume](double x, double y) { return volume->GetTopHeight(x, y); });
        terrain->DisallowCollisionsWith(shapes->GetCollisionFamily());
        terrain->Initialize(ChVector2d(0, 0));

        driver = chrono_types::make_shared<DelayedSpeedDriver>(rover_speed);
        viper = std::make_unique<Viper>(&sys, ViperWheelType::CylWheel);
        viper->SetDriver(driver);
        viper->SetWheelContactMaterial(contact_mat);
        const double x0 = work_site.min_x - 2.5;
        double ground = -1e9;
        for (double wx : {-0.64, 0.64})
            for (double wy : {-0.61, 0.61})
                ground = std::max(ground, height_at(x0 + wx, 1.7 + wy));
        viper->Initialize(ChFrame<>(ChVector3d(x0, 1.7, ground + wheel_radius + 0.05), QUNIT));
    }

    // The bucket, moved along the dig cycle
    const DigPlan plan = MakePlan(passes);
    DigCycle cycle(volume, plan);
    auto bucket = chrono_types::make_shared<ChBody>();
    bucket->SetFixed(true);
    bucket->EnableCollision(false);
    {
        const double t = 0.02;
        const ChVector3d s = bucket_size;
        auto steel = Material(0.75f, 0.55f, 0.15f);
        auto add = [&](const ChVector3d& size, const ChVector3d& at) {
            auto box = chrono_types::make_shared<ChVisualShapeBox>(size);
            box->AddMaterial(steel);
            bucket->AddVisualShape(box, ChFrame<>(at));
        };
        add(ChVector3d(s.x(), s.y(), t), ChVector3d(0, 0, -0.5 * s.z()));              // bottom
        add(ChVector3d(t, s.y(), s.z()), ChVector3d(-0.5 * s.x(), 0, 0));              // back
        add(ChVector3d(s.x(), t, s.z()), ChVector3d(0, -0.5 * s.y(), 0));              // sides
        add(ChVector3d(s.x(), t, s.z()), ChVector3d(0, 0.5 * s.y(), 0));
        add(ChVector3d(0.05, 0.05, 1.5), ChVector3d(-0.5 * s.x(), 0, 0.5 * s.z() + 0.75));  // stick
    }
    sys.AddBody(bucket);
    ChVector3d bucket_pos = cycle.Start().GetPos();
    double bucket_pitch = 0;
    bucket->SetPos(bucket_pos);

    ChExcavationTool tool(volume, bucket, chrono_types::make_shared<ChSdfBox>(bucket_size));
    ChExcavationTool::Soil soil;
    soil.gravity = gravity;
    soil.bulk_density = 1600;
    soil.cohesion = 1000;
    tool.SetSoil(soil);
    tool.SetCapacity(bucket_capacity);
    tool.SetRakeAngle(60 * CH_DEG_TO_RAD);

    // The soil the bucket carries and pours, as particles
    ChSoilParticles::Params pparams;
    pparams.spacing = particle_spacing;
    pparams.gravity = gravity;
    pparams.repose_angle = soil.repose_angle;
    auto particles = chrono_types::make_shared<ChSoilParticles>(&sys, volume, pparams);
    particles->SetUseVisualModel(!use_vis);
    particles->SetMaterial(Material(0.37f, 0.36f, 0.34f));
    const ChVector3d inner = 0.5 * bucket_size - ChVector3d(bucket_wall);
    tool.SetParticles(particles, ChAABB(-inner, inner));

    // Visualization
    SiteStats stats;
    stats.strokes = cycle.GetNumStrokes();
    stats.rover = use_rover;
    std::shared_ptr<vsg3d::ChVisualSystemVSG> vis;
    if (use_vis) {
        auto vis_planet = chrono_types::make_shared<ChPlanetVisualizationVSG>(world, site);
        vis_planet->AddHole(volume->GetHole());
        auto vis_site = chrono_types::make_shared<ChSiteVolumeVisualizationVSG>(shapes);
        auto vis_soil = chrono_types::make_shared<ChSoilParticlesVisualizationVSG>(particles);
        vis = chrono_types::make_shared<vsg3d::ChVisualSystemVSG>();
        vis->AttachSystem(&sys);
        vis->AttachPlugin(vis_planet);
        vis->AttachPlugin(vis_site);
        vis->AttachPlugin(vis_soil);
        vis->AddGuiComponent(chrono_types::make_shared<SiteStatsVSG>(stats));
        vis->SetWindowTitle("Excavation at a lunar work site");
        vis->SetWindowSize(1280, 800);
        vis->SetWindowPosition(100, 100);
        vis->SetBackgroundColor(ChColor(0.01f, 0.01f, 0.02f));
        if (close_camera)
            vis->AddCamera(ChVector3d(-3.4, -2.2, 1.5), ChVector3d(-1.6, 0.6, 0.1));
        else
            vis->AddCamera(ChVector3d(-3.5, -6.5, 3.5), ChVector3d(-0.6, 0.6, 0));
        vis->SetCameraVertical(CameraVerticalDir::Z);
        vis->SetCameraAngleDeg(45.0);
        vis->SetLightIntensity(1.0f);
        vis->SetLightDirection(1.2, 0.5);
        vis->Initialize();
    }

    const std::string out_dir = GetChronoOutputPath() + "PLANET_Excavation";
    if (snapshots && !CreateOutputDirectory(out_dir)) {
        std::cerr << "Error creating directory " << out_dir << std::endl;
        return 1;
    }
    double next_snapshot = 0;
    int snapshot = 0;

    // Simulation loop
    const auto wall_start = std::chrono::steady_clock::now();
    double mesh_time = 0;
    int step = 0;
    bool reported = false;
    while (sys.GetChTime() < end_time) {
        const double time = sys.GetChTime();
        if (vis && !vis->Run())
            break;

        // The bucket's next pose, and the soil it cuts on the way
        ChVector3d velocity;
        const bool dump = cycle.Advance(step_size, bucket_pos, bucket_pitch, velocity);
        bucket->SetPos(bucket_pos);
        bucket->SetRot(QuatFromAngleY(bucket_pitch));
        bucket->SetPosDt(velocity);
        tool.Update();
        if (dump)
            tool.Dump();  // what did not pour out falls
        particles->Advance(step_size);

        // The rover sets off when the digging is done
        if (use_rover && cycle.Done() && !driver->Started()) {
            driver->Start(time);
            std::cout << "Digging done at t = " << time << " s; the rover sets off" << std::endl;
        }

        if (step % render_steps == 0) {
            const auto t0 = std::chrono::steady_clock::now();
            stats.rebuilt = shapes->Update();
            volume->PublishTopSurface(*dug, site);
            particles->UpdateMesh();
            mesh_time += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

            stats.stroke = cycle.GetStroke();
            stats.depth = tool.GetCuttingDepth();
            stats.force = tool.GetResistanceForce().Length();
            stats.power = std::abs(tool.GetResistanceForce().Dot(velocity));
            stats.payload = tool.GetPayloadVolume();
            stats.carried = particles->GetNumCarried();
            stats.falling = particles->GetNumFalling();
            stats.resting = particles->GetNumResting();
            stats.excavated = tool.GetExcavatedVolume();
            stats.dumped = tool.GetDumpedVolume();
            stats.spilled = tool.GetSpilledVolume();
            stats.balance = volume->GetSoilVolume() + tool.GetPayloadVolume() + tool.GetPendingSpillVolume() + particles->GetFreeVolume() - initial_soil;
            stats.bricks = shapes->GetNumShapes();
            stats.triangles = shapes->GetNumTriangles();
            if (viper)
                stats.rover_x = viper->GetChassis()->GetPos().x();

            if (vis) {
                vis->BeginScene();
                vis->Render();
                vis->EndScene();
                if (snapshots && time >= next_snapshot) {
                    char name[64];
                    std::snprintf(name, sizeof(name), "/snapshot_%03d.png", snapshot++);
                    vis->WriteImageToFile(out_dir + name);
                    next_snapshot += snapshot_period;
                }
            }
        }

        sys.DoStepDynamics(step_size);
        if (viper) {
            viper->Update();
            terrain->UpdatePatch(viper->GetChassis()->GetPos());
        }

        if (++step % 500 == 0) {
            const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
            std::cout << "t = " << sys.GetChTime() << " s (" << sys.GetChTime() / wall << "x real time), stroke "
                      << std::min(cycle.GetStroke() + 1, stats.strokes) << "/" << stats.strokes << ", particles " << particles->GetNumCarried() << " carried "
                      << particles->GetNumFalling() << " falling, cut " << tool.GetExcavatedVolume() * 1e3
                      << " L, dumped " << tool.GetDumpedVolume() * 1e3 << " L, soil balance " << stats.balance * 1e3 << " L, "
                      << shapes->GetNumShapes() << " bricks, meshing " << mesh_time / wall * 100 << "% of the run";
            if (viper)
                std::cout << ", rover at (" << viper->GetChassis()->GetPos().x() << ", " << viper->GetChassis()->GetPos().y() << ", "
                          << viper->GetChassis()->GetPos().z() << ")";
            std::cout << std::endl;
        }
        if (cycle.Done() && !reported) {
            reported = true;
            std::cout << "Trench dug: " << tool.GetExcavatedVolume() * 1e3 << " L cut, " << tool.GetDumpedVolume() * 1e3 << " L dumped, "
                      << tool.GetSpilledVolume() * 1e3 << " L pushed ahead; soil in the site changed by "
                      << (volume->GetSoilVolume() - initial_soil) * 1e3 << " L; trench floor at "
                      << volume->GetTopHeight(0.0, 0.0) - volume->GetInitialHeight(0.0, 0.0) << " m" << std::endl;
            if (!use_rover && !vis)
                break;
        }
        if (use_rover && driver->Started() && viper->GetChassis()->GetPos().x() > work_site.max_x + 2.0 && !vis)
            break;
    }
    return 0;
}
