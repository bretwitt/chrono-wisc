#include "chrono_planet/shading/Hapke.h"

#include <algorithm>
#include <cmath>

#include "chrono_planet/core/MathUtil.h"

namespace qtplanet {
namespace {

constexpr float kPiF = 3.14159265358979f;

// Chandrasekhar H, Hapke's (2002) analytic approximation.
float chandrasekharH(float x, float w, float r0) {
    const float lx = std::log((1.0f + x) / std::max(x, 1e-4f));
    const float den = 1.0f - w * x * (r0 + (1.0f - 2.0f * r0 * x) * 0.5f * lx);
    return 1.0f / std::max(den, 1e-4f);
}

// Double Henyey-Greenstein. The backscatter lobe carries weight (1+c)/2.
float hapkePhase(float cosG, float b, float c) {
    const float b2 = b * b;
    const float back = (1.0f - b2) * std::pow(std::max(1.0f - 2.0f * b * cosG + b2, 1e-4f), -1.5f);
    const float fwd = (1.0f - b2) * std::pow(std::max(1.0f + 2.0f * b * cosG + b2, 1e-4f), -1.5f);
    return 0.5f * (1.0f + c) * back + 0.5f * (1.0f - c) * fwd;
}

// Hapke (1984) macroscopic roughness, effective cosines and the shadowing factor S.
struct RoughnessTerm {
    float mu0e, mue, S;
};

RoughnessTerm hapkeRoughness(float mu0, float mu, float cosPsi, float psi, float tanHalfG, float thetaBar) {
    const float t = std::tan(thetaBar);
    if (t < 1e-3f) {
        return {mu0, mu, 1.0f};   // smooth, the correction is identity
    }

    const float chi = 1.0f / std::sqrt(1.0f + kPiF * t * t);
    const float cotT = 1.0f / t;

    const float sinI = std::sqrt(std::max(1.0f - mu0 * mu0, 0.0f));
    const float sinE = std::sqrt(std::max(1.0f - mu * mu, 0.0f));
    const float cotI = mu0 / std::max(sinI, 1e-4f);
    const float cotE = mu / std::max(sinE, 1e-4f);

    const float E1i = std::exp(-2.0f / kPiF * cotT * cotI);
    const float E1e = std::exp(-2.0f / kPiF * cotT * cotE);
    const float E2i = std::exp(-1.0f / kPiF * cotT * cotT * cotI * cotI);
    const float E2e = std::exp(-1.0f / kPiF * cotT * cotT * cotE * cotE);

    const float sh = std::sin(psi * 0.5f);
    const float sin2 = sh * sh;                   // sin^2(psi/2)
    const float f = std::exp(-2.0f * tanHalfG);   // phase-angle decorrelation

    const float etaI = chi * (mu0 + sinI * t * E2i / std::max(2.0f - E1i, 1e-4f));
    const float etaE = chi * (mu + sinE * t * E2e / std::max(2.0f - E1e, 1e-4f));

    RoughnessTerm r;
    if (mu0 >= mu) {   // i <= e
        const float den = std::max(2.0f - E1e - (psi / kPiF) * E1i, 1e-4f);
        r.mu0e = chi * (mu0 + sinI * t * (cosPsi * E2e + sin2 * E2i) / den);
        r.mue = chi * (mu + sinE * t * (E2e - sin2 * E2i) / den);
        const float q = mu0 / std::max(etaI, 1e-4f);
        r.S = (r.mue / std::max(etaE, 1e-4f)) * q * chi / std::max(1.0f - f + f * chi * q, 1e-4f);
    } else {
        const float den = std::max(2.0f - E1i - (psi / kPiF) * E1e, 1e-4f);
        r.mu0e = chi * (mu0 + sinI * t * (E2i - sin2 * E2e) / den);
        r.mue = chi * (mu + sinE * t * (cosPsi * E2i + sin2 * E2e) / den);
        const float q = mu / std::max(etaE, 1e-4f);
        r.S = (r.mue / std::max(etaE, 1e-4f)) * (mu0 / std::max(etaI, 1e-4f)) * chi /
              std::max(1.0f - f + f * chi * q, 1e-4f);
    }
    return r;
}

}   // namespace

float thetaBarRadians(const HapkeParams& p) {
    return p.thetaBarDeg < 0.f ? -1.f : static_cast<float>(deg2rad(p.thetaBarDeg));
}

float porosityK(float phi) {
    const float x = 1.209f * std::pow(std::max(phi, 1e-4f), 2.0f / 3.0f);
    if (x >= 0.999f || x <= 1e-6f) {
        return 1.0f;   // outside the model's range
    }
    return -std::log(1.0f - x) / x;
}

float reflectance(float mu0, float mu, float cosG, float cosPsi, float psi, float thetaBar, const HapkeParams& p) {
    if (mu0 <= 0.0f) {
        return 0.0f;
    }
    const float tanHalfG = std::sqrt(std::max(1.0f - cosG, 0.0f) / std::max(1.0f + cosG, 1e-4f));
    const RoughnessTerm rg = hapkeRoughness(mu0, mu, cosPsi, psi, tanHalfG, thetaBar);

    const float K = porosityK(p.phi);
    const float gamma = std::sqrt(std::max(1.0f - p.w, 0.0f));
    const float r0 = (1.0f - gamma) / (1.0f + gamma);
    // Porosity enters the multiple-scattering term through the H arguments as well as the prefactor.
    const float M = chandrasekharH(rg.mu0e / K, p.w, r0) * chandrasekharH(rg.mue / K, p.w, r0) - 1.0f;

    const float Bs = 1.0f / (1.0f + tanHalfG / std::max(p.hs, 1e-4f));
    // Coherent backscatter (Hapke 2002 Eq. 32), multiplying the whole bracket. Bc(0) = 1 exactly.
    float Bc = 0.0f;
    if (p.Bc0 > 0.0f) {
        // Series below x ~ 1e-3 where (1 - e^-x)/x cancels in float, matching the GLSL mirror.
        const float x = tanHalfG / std::max(p.hc, 1e-4f);
        const float ex = (x < 1e-3f) ? (1.0f - 0.5f * x) : ((1.0f - std::exp(-x)) / x);
        Bc = (1.0f + ex) / (2.0f * (1.0f + x) * (1.0f + x));
    }
    const float ph = hapkePhase(cosG, p.b, p.c);
    const float LS = rg.mu0e / std::max(rg.mu0e + rg.mue, 1e-4f);
    return K * p.w * 0.25f / kPiF * LS * (ph * (1.0f + p.Bs0 * Bs) + M) * (1.0f + p.Bc0 * Bc) * rg.S;
}

float normalisation(const HapkeParams& p) {
    const float kDeg = kPiF / 180.0f;
    const float mu0 = std::cos(p.refIDeg * kDeg);
    const float mu = std::cos(p.refEDeg * kDeg);
    const float cosG = std::cos(p.refGDeg * kDeg);
    // At e = 0 the azimuth is undefined and irrelevant, and psi = 0 is the safe limit.
    const float r = reflectance(mu0, mu, cosG, 1.0f, 0.0f, p.nominalThetaBarDeg * kDeg, p);
    return (r > 1e-9f) ? r : 1e-9f;
}

}   // namespace qtplanet
