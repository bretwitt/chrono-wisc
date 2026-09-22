#ifndef QTPLANET_HAPKE_H
#define QTPLANET_HAPKE_H

#include "chrono_planet/ChApiPlanet.h"

// Hapke bidirectional reflectance for lunar regolith, after Batagoda et al. (arXiv:2410.04371, Eq. 9).
// Mirrored in shaders/common/hapke.glsl. The CPU copy evaluates the normalization once per frame.

namespace qtplanet {

// Model parameters. Table 1 of the paper, except w and Bs0, which are calibrated against measured lunar
// albedo and phase integral (vv/hapke groups 20-27). Do not restore the paper's values.
struct CH_PLANET_API HapkeParams {
    float w = 0.207f;      // single-scattering albedo, calibrated
    float b = 0.23955f;    // Henyey-Greenstein asymmetry
    float c = 0.30452f;    // backward/forward lobe weighting
    float Bs0 = 0.975f;    // shadow-hiding opposition amplitude, calibrated
    float hs = 0.07145f;   // shadow-hiding angular width, in tan(g/2)

    float Bc0 = 0.0f;       // coherent-backscatter amplitude, 0 reproduces the paper
    float hc = 0.017455f;   // coherent-backscatter angular width, in tan(g/2)

    float phi = 0.30f;   // filling factor, 1 - porosity

    // Photometric roughness. Negative derives it per pixel from the shader's slope variance.
    float thetaBarDeg = -1.0f;
    // Roughness used for the normalization only, which must be a scene constant.
    float nominalThetaBarDeg = 23.4f;

    // Geometry the color mosaic is photometrically normalized to.
    float refIDeg = 30.0f;
    float refEDeg = 0.0f;
    float refGDeg = 30.0f;
};

// thetaBarDeg in radians, preserving the negative sentinel.
float thetaBarRadians(const HapkeParams& p);

// Hapke (2008) porosity correction. phi = 0.3 gives K ~ 1.44.
float porosityK(float phi);

// Bidirectional reflectance r = f * mu0 in sr^-1, so L = r * J for collimated irradiance J.
// mu0 and mu are the incidence and emission cosines, cosG the phase cosine, psi the relative azimuth.
float reflectance(float mu0, float mu, float cosG, float cosPsi, float psi, float thetaBar,
                  const HapkeParams& p);

// Reflectance at the mosaic's reference geometry and nominal roughness. The shader divides by this.
float normalisation(const HapkeParams& p);

}   // namespace qtplanet

#endif   // QTPLANET_HAPKE_H
