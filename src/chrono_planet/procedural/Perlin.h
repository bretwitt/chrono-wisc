#ifndef QTPLANET_PERLIN_H
#define QTPLANET_PERLIN_H

#include "chrono_planet/ChApiPlanet.h"

// Classic 2D Perlin (1985) gradient noise with the quintic fade of Perlin (2002), and a sphere wrapper built from it.

#include <vector>

#include "chrono_planet/core/MathUtil.h"

namespace Perlin {

// Dimensionless noise, roughly [-1, 1]. x and y are lattice coordinates.
[[nodiscard]] float noise(float x, float y);

// Noise on a unit direction, the average of the three coordinate-plane
// samples. coordinateScale multiplies the unit coordinates into lattice units.
[[nodiscard]] float onSphere(const qtplanet::Vec3& unitDirection, float coordinateScale);

// Overwrites noiseValues[i] = noise(x[i], y[i]), bit-identical to scalar noise().
// sampleCount >= 0; each array has at least sampleCount elements and arrays do not overlap.
void writeNoiseRow(const float* x, const float* y, int sampleCount, float* noiseValues);

// Overwrites noiseValues[i] = onSphere(unitDirections[i], coordinateScale).
// sampleCount >= 0; input/output arrays have at least sampleCount elements.
// scratchValues is resized and overwritten; reuse across rows to avoid allocations.
// Input, output and scratch storage must not overlap.
void writeSphereNoiseRow(const qtplanet::Vec3* unitDirections, int sampleCount, float coordinateScale, float* noiseValues, std::vector<float>& scratchValues);

}   // namespace Perlin

#endif   // QTPLANET_PERLIN_H
