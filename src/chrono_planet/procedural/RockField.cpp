#include "chrono_planet/procedural/RockField.h"
#include "chrono_planet/core/SphereMath.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>

#include "chrono_planet/procedural/FieldHash.h"

using qtfield::fastFloor;
using qtfield::hash01;
using qtfield::kKmPerDeg;

namespace {

// Golombek and Rapp (1997) size-frequency distribution, F(D) = k exp(-q(k) D) with q(k) = kGrQa + kGrQb / k.
// k is the cumulative fractional area rocks cover.
constexpr double kRockCFA = 0.02;
constexpr double kGrQa = 1.79, kGrQb = 0.152;

// Class bin edges, meters of diameter, one octave per class.
constexpr double kBinEdgeM[] = {0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0};
constexpr int kOct = sizeof(kBinEdgeM) / sizeof(kBinEdgeM[0]) - 1;
static_assert(kOct == RockField::kClasses,
              "RockField::kClasses must match the bin table");

// Largest rock a class can hold, the bin's top edge.
constexpr double classDiamM(int octaveIndex) { return kBinEdgeM[octaveIndex + 1]; }

// Expected rocks per square meter in bin octaveIndex, Simpson's rule over the diameter density.
inline double binDensityPerM2(int octaveIndex) {
    const double exponentialDecay = kGrQa + kGrQb / kRockCFA;
    const double minDiameterM = kBinEdgeM[octaveIndex], maxDiameterM = kBinEdgeM[octaveIndex + 1];
    const auto densityAtDiameter = [exponentialDecay](double diameterM) {
        return (4.0 * kRockCFA * exponentialDecay / qtplanet::kPi) *
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
    double diameterDeg, cellSizeDeg, invCellSizeDeg;
    double density;                            // probability this cell holds a rock
    double invMinDiameterM, invMaxDiameterM;   // 1/D at the bin edges, for the within-bin draw
};
inline const Octave& octave(int octaveIndex) {
    static const auto octaveTable = []() {
        std::array<Octave, kOct> t{};
        for (int i = 0; i < kOct; ++i) {
            t[i].diameterDeg = (classDiamM(i) / 1000.0) / kKmPerDeg;
            t[i].cellSizeDeg = t[i].diameterDeg * kCellDiams;
            t[i].invCellSizeDeg = 1.0 / t[i].cellSizeDeg;
            // One candidate per cell, so the expected count must fit in it.
            const double cellM = classDiamM(i) * kCellDiams;
            t[i].density = std::min(1.0, binDensityPerM2(i) * cellM * cellM);
            t[i].invMinDiameterM = 1.0 / kBinEdgeM[i];
            t[i].invMaxDiameterM = 1.0 / kBinEdgeM[i + 1];
        }
        return t;
    }();
    return octaveTable[octaveIndex];
}

// Bedding octave fade, gated on the plinth's footprint.
inline double beddingFade(int octaveIndex, double invSpacingDeg) {
    return qtfield::resolutionFade(octave(octaveIndex).diameterDeg * kFadeDiams * invSpacingDeg, kFadeLoVerts, kFadeHiVerts);
}

// One rock's placement, a pure function of (octave, gx, gy). Returns nullopt if the cell holds no rock.
std::optional<RockInstance> makeInstance(int gx, int gy, int octaveIndex) {
    RockInstance rock;
    const Octave& octaveInfo = octave(octaveIndex);
    // Cells are square in degrees, so the per-square-meter density scales by 1/cos(lat).
    const double latDeg = (gy + 0.5) * octaveInfo.cellSizeDeg;
    const double cosLat = std::max(std::cos(qtplanet::deg2rad(latDeg)), kMinCellCosLat);
    if (hash01(gx, gy, octaveIndex * 29 + 101) > octaveInfo.density * cosLat) {
        return std::nullopt;
    }
    rock.id = (static_cast<std::uint64_t>(octaveIndex) << 60) ^
              (static_cast<std::uint64_t>(static_cast<std::uint32_t>(gx)) << 30) ^
              static_cast<std::uint64_t>(static_cast<std::uint32_t>(gy));
    rock.lonDeg = (gx + hash01(gx, gy, octaveIndex * 29 + 102)) * octaveInfo.cellSizeDeg;
    rock.latDeg = (gy + hash01(gx, gy, octaveIndex * 29 + 103)) * octaveInfo.cellSizeDeg;
    // Diameter within the bin from the D^-2 law, continuous across bin edges.
    const double u = hash01(gx, gy, octaveIndex * 29 + 104);
    const double diameterM = 1.0 / (octaveInfo.invMinDiameterM + u * (octaveInfo.invMaxDiameterM - octaveInfo.invMinDiameterM));
    rock.radiusM = static_cast<float>(0.5 * diameterM);
    rock.sizeClass = static_cast<std::uint8_t>(octaveIndex);

    rock.yawRad = static_cast<float>(qtplanet::kTwoPi * hash01(gx, gy, octaveIndex * 29 + 110));
    const double tiltSample = hash01(gx, gy, octaveIndex * 29 + 111);
    rock.tiltRad = static_cast<float>(kMaxTiltRad * tiltSample * tiltSample);
    rock.tiltAzRad = static_cast<float>(qtplanet::kTwoPi * hash01(gx, gy, octaveIndex * 29 + 112));
    rock.buryFrac = static_cast<float>(0.15 + 0.25 * hash01(gx, gy, octaveIndex * 29 + 113));
    // Resolved modulo RockMeshes::kCount by the consumer.
    rock.meshId = static_cast<std::uint16_t>(hash01(gx, gy, octaveIndex * 29 + 114) * 65535.f);
    return rock;
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

}   // namespace

double RockField::classDiameterM(int sizeClass) {
    return (sizeClass >= 0 && sizeClass < kOct) ? classDiamM(sizeClass) : 0.0;
}

// Coarsest class, not class 0, since the bin table runs finest first.
double RockField::maxRadiusM() { return 0.5 * classDiamM(kOct - 1); }

std::vector<RockInstance> RockField::query(double minLonDeg, double minLatDeg, double maxLonDeg, double maxLatDeg,
                                           double minRadiusM) {
    std::vector<RockInstance> out;
    for (int octaveIndex = 0; octaveIndex < kOct; ++octaveIndex) {
        // Whole class too small to be asked for.
        if (0.5 * classDiamM(octaveIndex) < minRadiusM) {
            continue;
        }
        const Octave& octaveInfo = octave(octaveIndex);
        const int gx0 = fastFloor(minLonDeg * octaveInfo.invCellSizeDeg), gx1 = fastFloor(maxLonDeg * octaveInfo.invCellSizeDeg);
        const int gy0 = fastFloor(minLatDeg * octaveInfo.invCellSizeDeg), gy1 = fastFloor(maxLatDeg * octaveInfo.invCellSizeDeg);
        for (int gy = gy0; gy <= gy1; ++gy) {
            for (int gx = gx0; gx <= gx1; ++gx) {
                const auto instance = makeInstance(gx, gy, octaveIndex);
                if (!instance) {
                    continue;
                }
                const RockInstance& rock = *instance;
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

std::array<double, 4> RockField::orientation(const RockInstance& rock) {
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
    return qtplanet::normalized4(orientation);
}

double RockField::centreRiseM(const RockInstance& rock, double bottomExtentM,
                              double totalHeightM) {
    // Bury at least kMinBuryM, and never past the equator.
    const double burialDepthM = std::min(0.5 * totalHeightM,
                                         std::max(kMinBuryM, rock.buryFrac * totalHeightM));
    return bottomExtentM - burialDepthM;
}

// One point. Only the neighboring cell on the side the point sits toward can contribute.
double RockField::beddingAt(double lonDeg, double latDeg, double sampleSpacingDeg) {
    const double invSampleSpacingDeg = 1.0 / std::max(sampleSpacingDeg, 1e-9);
    const double cosLat = qtplanet::cosLatClamped(latDeg);
    const double mPerDeg = kKmPerDeg * 1000.0;
    double reliefM = 0.0;
    for (int octaveIndex = 0; octaveIndex < kOct; ++octaveIndex) {
        const double fade = beddingFade(octaveIndex, invSampleSpacingDeg);
        if (fade <= 0.0) {
            break;
        }
        const Octave& octaveInfo = octave(octaveIndex);
        const double cx = lonDeg * octaveInfo.invCellSizeDeg, cy = latDeg * octaveInfo.invCellSizeDeg;
        const int ix = fastFloor(cx), iy = fastFloor(cy);
        const double fx = cx - ix, fy = cy - iy;
        const int ox0 = fx < kSkipCells ? -1 : 0, ox1 = fx > 1.0 - kSkipCells ? 1 : 0;
        const int oy0 = fy < kSkipCells ? -1 : 0, oy1 = fy > 1.0 - kSkipCells ? 1 : 0;
        for (int oy = oy0; oy <= oy1; ++oy) {
            for (int ox = ox0; ox <= ox1; ++ox) {
                const auto instance = makeInstance(ix + ox, iy + oy, octaveIndex);
                if (!instance) {
                    continue;
                }
                const RockInstance& rock = *instance;
                double dLon = lonDeg - rock.lonDeg;
                if (dLon > 180.0) {
                    dLon -= 360.0;
                } else if (dLon < -180.0) {
                    dLon += 360.0;
                }
                const double dx = dLon * cosLat * mPerDeg, dy = (latDeg - rock.latDeg) * mPerDeg;
                reliefM += Bedding(rock.radiusM).reliefMAtDistanceSquared(dx * dx + dy * dy) * fade;
            }
        }
    }
    return reliefM;
}

// Rows [rowBegin, rowEnd) only. Walks every cell whose plinth can reach the band, so sums match serial bit for bit.
void RockField::addBeddingGridRows(const qtfield::GridSpec& grid, qtfield::RowRange rows, std::vector<double>& heightsM) {
    const auto [originLonDeg, originLatDeg, stepLonDeg, stepLatDeg, sampleCount, sampleSpacingDeg] = grid;
    const auto [rowBegin, rowEnd] = rows;
    if (rowBegin >= rowEnd) {
        return;
    }
    const double invSampleSpacingDeg = 1.0 / std::max(sampleSpacingDeg, 1e-9);
    const double lon1 = originLonDeg + (sampleCount - 1) * stepLonDeg;
    const double latB0 = originLatDeg + rowBegin * stepLatDeg, latB1 = originLatDeg + (rowEnd - 1) * stepLatDeg;
    const double mPerDeg = kKmPerDeg * 1000.0;
    std::vector<double> cosLat(sampleCount);
    for (int j = rowBegin; j < rowEnd; ++j) {
        cosLat[j] = qtplanet::cosLatClamped(originLatDeg + j * stepLatDeg);
    }
    const double invStepLon = 1.0 / stepLonDeg, invStepLat = 1.0 / stepLatDeg;

    for (int octaveIndex = 0; octaveIndex < kOct; ++octaveIndex) {
        const double fade = beddingFade(octaveIndex, invSampleSpacingDeg);
        if (fade <= 0.0) {
            break;
        }
        const Octave& octaveInfo = octave(octaveIndex);
        // Cells whose plinth can touch the band.
        const int gx0 = fastFloor((originLonDeg - kReachCells * octaveInfo.cellSizeDeg) * octaveInfo.invCellSizeDeg);
        const int gx1 = fastFloor((lon1 + kReachCells * octaveInfo.cellSizeDeg) * octaveInfo.invCellSizeDeg);
        const int gy0 = fastFloor((latB0 - kReachCells * octaveInfo.cellSizeDeg) * octaveInfo.invCellSizeDeg);
        const int gy1 = fastFloor((latB1 + kReachCells * octaveInfo.cellSizeDeg) * octaveInfo.invCellSizeDeg);
        for (int gy = gy0; gy <= gy1; ++gy) {
            for (int gx = gx0; gx <= gx1; ++gx) {
                const auto instance = makeInstance(gx, gy, octaveIndex);
                if (!instance) {
                    continue;
                }
                const RockInstance& rock = *instance;
                const double reachM = kBeddingRadiusFrac * rock.radiusM;
                const double reachLat = reachM / mPerDeg;
                const Bedding bed(rock.radiusM);
                const int j0 = std::max(rowBegin, static_cast<int>(std::ceil((rock.latDeg - reachLat - originLatDeg) * invStepLat)));
                const int j1 = std::min(rowEnd - 1, static_cast<int>(std::floor((rock.latDeg + reachLat - originLatDeg) * invStepLat)));
                for (int j = j0; j <= j1; ++j) {
                    const double dy = (originLatDeg + j * stepLatDeg - rock.latDeg) * mPerDeg;
                    const double reachLon = reachM / (mPerDeg * cosLat[j]);
                    const int i0 = std::max(0, static_cast<int>(std::ceil((rock.lonDeg - reachLon - originLonDeg) * invStepLon)));
                    const int i1 = std::min(sampleCount - 1, static_cast<int>(std::floor((rock.lonDeg + reachLon - originLonDeg) * invStepLon)));
                    double* row = &heightsM[static_cast<size_t>(j) * sampleCount];
                    const double rowCosLat = cosLat[j];
                    for (int i = i0; i <= i1; ++i) {
                        const double dx = (originLonDeg + i * stepLonDeg - rock.lonDeg) * rowCosLat * mPerDeg;
                        row[i] += bed.reliefMAtDistanceSquared(dx * dx + dy * dy) * fade;
                    }
                }
            }
        }
    }
}
