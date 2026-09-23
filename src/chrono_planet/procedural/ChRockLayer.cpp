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

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <stdexcept>

#include "chrono_planet/core/Parallel.h"
#include "chrono_planet/core/SphereMath.h"
#include "chrono_planet/procedural/ChRockLayer.h"
#include "chrono_planet/procedural/FieldGrid.h"
#include "chrono_planet/core/FieldHash.h"

namespace chrono {
namespace planet {

using field::fastFloor;
using field::hash01;
using field::seeded;

namespace {

// Bedding, the shallow plinth a rock sits in, a capped fraction of the radius tall and kBeddingRadiusFrac wider.
constexpr double kBeddingRiseFrac = 0.12;
constexpr double kBeddingRadiusFrac = 1.45;
constexpr double kBeddingMaxRiseM = 0.30;
// Jitter-grid cell size in diameters of the class's largest rock.
constexpr double kCellDiams = 3.0;
// A plinth's farthest reach in cells, a bound on kBeddingRadiusFrac over twice kCellDiams.
constexpr double kReachCells = 0.25;
// Neighbor cells whose near edge is farther than this from a point are skipped, generous over kReachCells.
constexpr double kSkipCells = 0.4;
// Plinth footprint in rock diameters for the octave fade, and the fade bounds in vertices across it.
constexpr double kFadeDiams = 1.2;
constexpr double kFadeLoVerts = 2.0, kFadeHiVerts = 4.0;
// Least burial depth (meters), so the mesh dipping between vertices cannot show daylight under an edge.
constexpr double kMinBuryM = 0.10;
// Lean off vertical, squared to skew toward upright.
constexpr double kMaxTiltRad = 0.60;
// Floor on cos(lat) for the per-cell density, looser than the grid walkers' so polar cells still thin out.
constexpr double kMinCellCosLat = 0.02;

struct Octave {
    double minDiameterM, maxDiameterM;
    double diameterDeg, cellSizeDeg, invCellSizeDeg;
    double density;                            // probability this cell holds a rock
    double invMinDiameterM, invMaxDiameterM;   // 1/D at the bin edges, for the within-bin draw
};

// Expected rocks per square meter with diameters in [minDiameterM, maxDiameterM], Simpson's rule over the diameter density.
double binDensityPerM2(double coverage, double qa, double qb, double minDiameterM, double maxDiameterM) {
    const double exponentialDecay = qa + qb / coverage;
    const auto densityAtDiameter = [&](double diameterM) {
        return (4.0 * coverage * exponentialDecay / util::kPi) *
               std::exp(-exponentialDecay * diameterM) / (diameterM * diameterM);
    };
    constexpr int kSteps = 64;   // even
    const double diameterStepM = (maxDiameterM - minDiameterM) / kSteps;
    double sum = densityAtDiameter(minDiameterM) + densityAtDiameter(maxDiameterM);
    for (int i = 1; i < kSteps; ++i) {
        sum += densityAtDiameter(minDiameterM + i * diameterStepM) * ((i & 1) ? 4.0 : 2.0);
    }
    return sum * diameterStepM / 3.0;
}

// A rock's plinth, squared reach and rise, hoisted out of the per-sample evaluation.
struct Bedding {
    double radiusSquaredM2, riseM;
    explicit Bedding(double radiusM) {
        const double radiusMWithBedding = kBeddingRadiusFrac * radiusM;
        radiusSquaredM2 = radiusMWithBedding * radiusMWithBedding;
        riseM = std::min(kBeddingMaxRiseM, kBeddingRiseFrac * radiusM);
    }
    // Relief in meters at a squared distance in square meters from the center. Straight-line so a footprint row vectorizes.
    double reliefMAtDistanceSquared(double distanceSquaredM2) const {
        const double u = 1.0 - distanceSquaredM2 / radiusSquaredM2;
        return (distanceSquaredM2 < radiusSquaredM2) ? riseM * u * u : 0.0;   // C1, zero slope at the rim
    }
};

}  // namespace

// Parameters with their derived per-class tables.
struct ChRockLayer::Model {
    Params params;
    double kmPerDeg;
    std::vector<Octave> octaves;

    int count() const { return static_cast<int>(octaves.size()); }
    float hash(int gx, int gy, int salt) const { return hash01(gx, gy, seeded(salt, params.seed)); }

    // Bedding octave fade, gated on the plinth's footprint.
    double beddingFade(int octaveIndex, double invSpacingDeg) const {
        return field::resolutionFade(octaves[octaveIndex].diameterDeg * kFadeDiams * invSpacingDeg, kFadeLoVerts, kFadeHiVerts);
    }

    // One rock's placement, a pure function of (octave, gx, gy). Returns nullopt if the cell holds no rock.
    std::optional<ChRockInstance> makeInstance(int gx, int gy, int octaveIndex) const {
        ChRockInstance rock;
        const Octave& octaveInfo = octaves[octaveIndex];
        // Cells are square in degrees, so the per-square-meter density scales by 1/cos(lat).
        const double latDeg = (gy + 0.5) * octaveInfo.cellSizeDeg;
        const double cosLat = std::max(std::cos(util::deg2rad(latDeg)), kMinCellCosLat);
        if (hash(gx, gy, octaveIndex * 29 + 101) > octaveInfo.density * cosLat) {
            return std::nullopt;
        }
        rock.id = (static_cast<std::uint64_t>(octaveIndex) << 60) ^
                  (static_cast<std::uint64_t>(static_cast<std::uint32_t>(gx)) << 30) ^
                  static_cast<std::uint64_t>(static_cast<std::uint32_t>(gy));
        rock.lonDeg = (gx + hash(gx, gy, octaveIndex * 29 + 102)) * octaveInfo.cellSizeDeg;
        rock.latDeg = (gy + hash(gx, gy, octaveIndex * 29 + 103)) * octaveInfo.cellSizeDeg;
        // Diameter within the bin from the D^-2 law, continuous across bin edges.
        const double u = hash(gx, gy, octaveIndex * 29 + 104);
        const double diameterM = 1.0 / (octaveInfo.invMinDiameterM + u * (octaveInfo.invMaxDiameterM - octaveInfo.invMinDiameterM));
        rock.radiusM = static_cast<float>(0.5 * diameterM);
        rock.sizeClass = static_cast<std::uint8_t>(octaveIndex);

        rock.yawRad = static_cast<float>(util::kTwoPi * hash(gx, gy, octaveIndex * 29 + 110));
        const double tiltSample = hash(gx, gy, octaveIndex * 29 + 111);
        rock.tiltRad = static_cast<float>(kMaxTiltRad * tiltSample * tiltSample);
        rock.tiltAzRad = static_cast<float>(util::kTwoPi * hash(gx, gy, octaveIndex * 29 + 112));
        rock.buryFrac = static_cast<float>(0.15 + 0.25 * hash(gx, gy, octaveIndex * 29 + 113));
        // Resolved modulo RockMeshes::kCount by the consumer.
        rock.meshId = static_cast<std::uint16_t>(hash(gx, gy, octaveIndex * 29 + 114) * 65535.f);
        return rock;
    }
};

// -----------------------------------------------------------------------------

ChRockLayer::ChRockLayer(const ChPlanetBody& body, const Params& params) {
    const auto& edges = params.bin_edges_m;
    for (size_t i = 0; i < edges.size(); ++i) {
        if (!(edges[i] > 0) || (i > 0 && !(edges[i] > edges[i - 1])))
            throw std::invalid_argument("ChRockLayer: bin_edges_m must be positive and strictly increasing");
    }
    if (params.coverage < 0 || params.coverage >= 1)
        throw std::invalid_argument("ChRockLayer: coverage must lie in [0, 1)");

    auto model = std::make_unique<Model>();
    model->params = params;
    model->kmPerDeg = field::kmPerDeg(body.GetRadius());
    // With no coverage there are no rocks at all, so no classes to walk.
    if (params.coverage > 0) {
        for (size_t i = 0; i + 1 < edges.size(); ++i) {
            Octave o;
            o.minDiameterM = edges[i];
            o.maxDiameterM = edges[i + 1];
            o.diameterDeg = (o.maxDiameterM / 1000.0) / model->kmPerDeg;
            o.cellSizeDeg = o.diameterDeg * kCellDiams;
            o.invCellSizeDeg = 1.0 / o.cellSizeDeg;
            // One candidate per cell, so the expected count must fit in it.
            const double cellM = o.maxDiameterM * kCellDiams;
            o.density = std::min(1.0, binDensityPerM2(params.coverage, params.qa, params.qb, o.minDiameterM, o.maxDiameterM) * cellM * cellM);
            o.invMinDiameterM = 1.0 / o.minDiameterM;
            o.invMaxDiameterM = 1.0 / o.maxDiameterM;
            model->octaves.push_back(o);
        }
    }
    m_model = std::move(model);
}

ChRockLayer::~ChRockLayer() = default;

const ChRockLayer::Params& ChRockLayer::GetParams() const {
    return m_model->params;
}

int ChRockLayer::GetNumClasses() const {
    return m_model->count();
}

double ChRockLayer::GetClassDiameter(int sizeClass) const {
    return (sizeClass >= 0 && sizeClass < m_model->count()) ? m_model->octaves[sizeClass].maxDiameterM : 0.0;
}

// Coarsest class, not class 0, since the bins run finest first.
double ChRockLayer::GetMaxRadius() const {
    return m_model->octaves.empty() ? 0.0 : 0.5 * m_model->octaves.back().maxDiameterM;
}

std::vector<ChRockInstance> ChRockLayer::Query(double minLonDeg, double minLatDeg, double maxLonDeg, double maxLatDeg,
                                               double minRadiusM) const {
    const Model& m = *m_model;
    std::vector<ChRockInstance> out;
    for (int octaveIndex = 0; octaveIndex < m.count(); ++octaveIndex) {
        const Octave& octaveInfo = m.octaves[octaveIndex];
        // Whole class too small to be asked for.
        if (0.5 * octaveInfo.maxDiameterM < minRadiusM) {
            continue;
        }
        const int gx0 = fastFloor(minLonDeg * octaveInfo.invCellSizeDeg), gx1 = fastFloor(maxLonDeg * octaveInfo.invCellSizeDeg);
        const int gy0 = fastFloor(minLatDeg * octaveInfo.invCellSizeDeg), gy1 = fastFloor(maxLatDeg * octaveInfo.invCellSizeDeg);
        for (int gy = gy0; gy <= gy1; ++gy) {
            for (int gx = gx0; gx <= gx1; ++gx) {
                const auto instance = m.makeInstance(gx, gy, octaveIndex);
                if (!instance) {
                    continue;
                }
                const ChRockInstance& rock = *instance;
                if (rock.radiusM < minRadiusM) {
                    continue;
                }
                // Half-open on the max edges, so adjacent queries neither duplicate nor drop a rock.
                if (rock.lonDeg < minLonDeg || rock.lonDeg >= maxLonDeg) {
                    continue;
                }
                if (rock.latDeg < minLatDeg || rock.latDeg >= maxLatDeg) {
                    continue;
                }
                out.push_back(rock);
            }
        }
    }
    return out;
}

std::array<double, 4> ChRockLayer::GetOrientation(const ChRockInstance& rock) {
    // Yaw about the local up, then a lean about a horizontal axis. Double trig keeps the quaternion unit.
    const double halfYawRad = 0.5 * static_cast<double>(rock.yawRad);
    const double halfTiltRad = 0.5 * static_cast<double>(rock.tiltRad);
    const double tiltAzimuthRad = static_cast<double>(rock.tiltAzRad);
    const double sinHalfYaw = std::sin(halfYawRad), cosHalfYaw = std::cos(halfYawRad);
    const double sinHalfTilt = std::sin(halfTiltRad), cosHalfTilt = std::cos(halfTiltRad);
    // The tilt axis is horizontal, perpendicular to the lean azimuth.
    const double tiltAxisX = -std::sin(tiltAzimuthRad), tiltAxisY = std::cos(tiltAzimuthRad);
    const double tiltQuaternion[4] = {tiltAxisX * sinHalfTilt, tiltAxisY * sinHalfTilt, 0.0, cosHalfTilt};   // (x, y, z, w)
    const double yawQuaternion[4] = {0.0, 0.0, sinHalfYaw, cosHalfYaw};
    // Hamilton product tiltQuaternion (x) yawQuaternion.
    const std::array<double, 4> orientation{
        tiltQuaternion[3] * yawQuaternion[0] + tiltQuaternion[0] * yawQuaternion[3] + tiltQuaternion[1] * yawQuaternion[2] - tiltQuaternion[2] * yawQuaternion[1],
        tiltQuaternion[3] * yawQuaternion[1] - tiltQuaternion[0] * yawQuaternion[2] + tiltQuaternion[1] * yawQuaternion[3] + tiltQuaternion[2] * yawQuaternion[0],
        tiltQuaternion[3] * yawQuaternion[2] + tiltQuaternion[0] * yawQuaternion[1] - tiltQuaternion[1] * yawQuaternion[0] + tiltQuaternion[2] * yawQuaternion[3],
        tiltQuaternion[3] * yawQuaternion[3] - tiltQuaternion[0] * yawQuaternion[0] - tiltQuaternion[1] * yawQuaternion[1] - tiltQuaternion[2] * yawQuaternion[2]};
    return util::normalized4(orientation);
}

double ChRockLayer::GetCenterRise(const ChRockInstance& rock, double bottomExtentM, double totalHeightM) {
    // Bury at least kMinBuryM, and never past the equator.
    const double burialDepthM = std::min(0.5 * totalHeightM, std::max(kMinBuryM, rock.buryFrac * totalHeightM));
    return bottomExtentM - burialDepthM;
}

// One point. Only the neighboring cell on the side the point sits toward can contribute.
double ChRockLayer::GetHeight(double lonDeg, double latDeg, double sampleSpacingDeg) const {
    const Model& m = *m_model;
    const double invSampleSpacingDeg = 1.0 / std::max(sampleSpacingDeg, 1e-9);
    const double cosLat = util::cosLatClamped(latDeg);
    const double mPerDeg = m.kmPerDeg * 1000.0;
    double reliefM = 0.0;
    for (int octaveIndex = 0; octaveIndex < m.count(); ++octaveIndex) {
        const double fade = m.beddingFade(octaveIndex, invSampleSpacingDeg);
        if (fade <= 0.0) {
            break;
        }
        const Octave& octaveInfo = m.octaves[octaveIndex];
        field::forEachCellNear(lonDeg * octaveInfo.invCellSizeDeg, latDeg * octaveInfo.invCellSizeDeg, kSkipCells, [&](int gx, int gy) {
            const auto instance = m.makeInstance(gx, gy, octaveIndex);
            if (!instance) {
                return;
            }
            const ChRockInstance& rock = *instance;
            double dLon = lonDeg - rock.lonDeg;
            if (dLon > 180.0) {
                dLon -= 360.0;
            } else if (dLon < -180.0) {
                dLon += 360.0;
            }
            const double dx = dLon * cosLat * mPerDeg, dy = (latDeg - rock.latDeg) * mPerDeg;
            reliefM += Bedding(rock.radiusM).reliefMAtDistanceSquared(dx * dx + dy * dy) * fade;
        });
    }
    return reliefM;
}

namespace {

// Rows `rows` only. Walks every cell whose plinth can reach the band, so sums match serial bit for bit.
void addBeddingGridRows(const ChRockLayer::Model& m, const ChGeoGrid& grid, field::RowRange rows, std::vector<double>& heightsM) {
    if (rows.begin >= rows.end) {
        return;
    }
    const double invSampleSpacingDeg = 1.0 / std::max(grid.spacing, 1e-9);
    const double lon1 = grid.lon0 + (grid.n - 1) * grid.step_lon;
    const double latB0 = grid.lat0 + rows.begin * grid.step_lat, latB1 = grid.lat0 + (rows.end - 1) * grid.step_lat;
    const double mPerDeg = m.kmPerDeg * 1000.0;
    std::vector<double> cosLat(grid.n);
    for (int j = rows.begin; j < rows.end; ++j) {
        cosLat[j] = util::cosLatClamped(grid.lat0 + j * grid.step_lat);
    }

    for (int octaveIndex = 0; octaveIndex < m.count(); ++octaveIndex) {
        const double fade = m.beddingFade(octaveIndex, invSampleSpacingDeg);
        if (fade <= 0.0) {
            break;
        }
        const Octave& octaveInfo = m.octaves[octaveIndex];
        // Cells whose plinth can touch the band.
        const int gx0 = fastFloor((grid.lon0 - kReachCells * octaveInfo.cellSizeDeg) * octaveInfo.invCellSizeDeg);
        const int gx1 = fastFloor((lon1 + kReachCells * octaveInfo.cellSizeDeg) * octaveInfo.invCellSizeDeg);
        const int gy0 = fastFloor((latB0 - kReachCells * octaveInfo.cellSizeDeg) * octaveInfo.invCellSizeDeg);
        const int gy1 = fastFloor((latB1 + kReachCells * octaveInfo.cellSizeDeg) * octaveInfo.invCellSizeDeg);
        for (int gy = gy0; gy <= gy1; ++gy) {
            for (int gx = gx0; gx <= gx1; ++gx) {
                const auto instance = m.makeInstance(gx, gy, octaveIndex);
                if (!instance) {
                    continue;
                }
                const ChRockInstance& rock = *instance;
                const double reachM = kBeddingRadiusFrac * rock.radiusM;
                const Bedding bed(rock.radiusM);
                field::scatterFootprint(grid, rows, rock.lonDeg, rock.latDeg, reachM / mPerDeg,
                                        [&](int j) { return reachM / (mPerDeg * cosLat[j]); },
                                        [&](int j, int i0, int i1) {
                                            const double dy = (grid.lat0 + j * grid.step_lat - rock.latDeg) * mPerDeg;
                                            // Locals, since through `grid` they could alias the row and the loop stops vectorizing.
                                            const double rowCosLat = cosLat[j], lon0 = grid.lon0, stepLon = grid.step_lon, rockLon = rock.lonDeg;
                                            double* row = &heightsM[static_cast<size_t>(j) * grid.n];
                                            for (int i = i0; i <= i1; ++i) {
                                                const double dx = (lon0 + i * stepLon - rockLon) * rowCosLat * mPerDeg;
                                                row[i] += bed.reliefMAtDistanceSquared(dx * dx + dy * dy) * fade;
                                            }
                                        });
            }
        }
    }
}

}  // namespace

void ChRockLayer::AddToGrid(const ChGeoGrid& grid, std::vector<double>& heights) const {
    const Model& m = *m_model;
    if (m.octaves.empty()) {
        return;
    }
    field::forEachRowBand(grid.n, [&](field::RowRange rows) { addBeddingGridRows(m, grid, rows, heights); });
}

}  // namespace planet
}  // namespace chrono
