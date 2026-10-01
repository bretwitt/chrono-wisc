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
// VIPER rover on hybrid soil (vehicle::PlanetHybridTerrain) at a lunar landing
// site, seen from the side by a Chrono::Sensor camera that follows it. The
// rover drives down into a shallow valley on SCM soil, which runs about as fast
// as real time; when the ground ahead rises past 10 degrees, a CRM window is
// seeded around it from the ruts SCM left, warms up without load, and takes
// the wheels' load over 0.3 s. On the slope the CRM soil gives the wheels far
// less traction than SCM did: slip rises from about 0.1 to 0.5 and more, and
// the rover slows to a crawl, as in demo_PLANET_Viper_CRM_HillClimb.
//
// With --save --no-window, the camera's frames (10 per second) go to
// DEMO_OUTPUT/PLANET_Viper_Hybrid_HillClimb/side, and soil_mode.txt there
// gives the soil model at each frame (time, then mode|CRM share|slope ahead
// (deg)|slip). --soil scm keeps SCM throughout and --soil crm uses CRM
// throughout, for comparison. See PlanetViperCRMSensor.h for the other flags.
//
// =============================================================================

#include "PlanetViperCRMSensor.h"

int main(int argc, char* argv[]) {
    ViperCRMScenario scenario;
    scenario.name = "PLANET_Viper_Hybrid_HillClimb";
    scenario.start_x = 5.45, scenario.start_y = -5.8, scenario.start_heading = 285;  // down into the valley, facing the hill
    scenario.speed = 0.8;
    scenario.end_time = 50;
    scenario.rate = 10;
    scenario.camera = "side";
    scenario.sun_elevation = 30;
    scenario.sun_azimuth = 165;
    scenario.soil = "hybrid";
    return RunViperCRMSensor(argc, argv, scenario);
}
