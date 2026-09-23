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
// Procedural impact craters with a configurable size-frequency distribution
// and depth/rim morphometry.
//
// =============================================================================

#ifndef CH_CRATER_LAYER_H
#define CH_CRATER_LAYER_H

#include <memory>
#include <vector>

#include "chrono_planet/ChPlanetBody.h"
#include "chrono_planet/procedural/ChReliefLayer.h"

namespace chrono {
namespace planet {

/// @addtogroup planet_module
/// @{

/// Impact craters scattered over the body, a pure function of longitude and latitude.
/// Craters come in size classes. Each class has a jittered grid whose cells are 2.5 class diameters
/// wide, and a cell holds a crater with the class's density (probability). A crater's diameter
/// varies within its class and its age varies from fresh (sharp rim, ejecta rays) to subdued.
/// Depth and rim height follow power laws in diameter, with a simple-to-complex transition above
/// which floors flatten and central peaks appear (the form of Pike, 1977).
///
/// The default parameters describe no craters. See planets/moon for a calibrated lunar set, and use
/// it as a template for other bodies: the transition diameter scales roughly with 1/gravity.
class CH_PLANET_API ChCraterLayer : public ChReliefLayer {
  public:
    struct Params {
        /// Class diameters (km), strictly decreasing.
        std::vector<double> diameters_km;
        /// Per-class probability, in [0, 1], that a jitter cell holds a crater. Same length as diameters_km.
        std::vector<double> densities;

        double complex_diameter_km = 0;  ///< simple-to-complex transition diameter (km)
        double simple_depth_ratio = 0;   ///< simple crater depth / diameter
        double simple_rim_ratio = 0;     ///< simple crater rim height / diameter
        double complex_depth_coef = 0;   ///< complex crater depth = coef * D^exp (km)
        double complex_depth_exp = 0;    ///< exponent of the complex depth law
        double complex_rim_coef = 0;     ///< complex crater rim height = coef * D^exp (km)
        double complex_rim_exp = 0;      ///< exponent of the complex rim law
        double floor_flattening_km = 30; ///< diameter range above the transition over which floors flatten fully (km)

        int seed = 0;  ///< varies the realization; 0 reproduces the reference field
    };

    /// Construct the layer for a body. Throws std::invalid_argument on inconsistent parameters.
    ChCraterLayer(const ChPlanetBody& body, const Params& params);
    ~ChCraterLayer();

    const Params& GetParams() const;

    virtual double GetHeight(double lon_deg, double lat_deg, double spacing_deg) const override;
    virtual void AddToGrid(const ChGeoGrid& grid, std::vector<double>& heights) const override;

    struct Model;

  private:
    std::unique_ptr<const Model> m_model;
};

/// @} planet_module

}  // namespace planet
}  // namespace chrono

#endif
