#include "chrono_planet/procedural/CraterField.h"
#include "chrono_planet/core/SphereMath.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <utility>
#include <vector>

#include "chrono_planet/core/Parallel.h"
#include "chrono_planet/procedural/FieldHash.h"
#include "chrono_planet/procedural/RockField.h"
#include "chrono_planet/procedural/SurfaceRoughness.h"

using qtfield::bump;
using qtfield::fastFloor;
using qtfield::hash01;
using qtfield::kKmPerDeg;

namespace {

// Rim bump support either side of the crest, in rim sigmas.
constexpr double kRimWidth = 2.2;
// Rim half-width cap in radiusFraction. Past kReachQ the profile is identically zero.
constexpr double kRimSigmaMax = 0.3;
constexpr double kReachQ = 1.0 + kRimWidth * kRimSigmaMax;
// Ejecta apron edge (in radiusFraction), also a crater's farthest reach.
constexpr double kEjectaQ = 2.0;
// Central peak diameter and height as fractions of the crater's own, after Pike.
constexpr double kPeakDiamFrac = 0.22;
constexpr double kPeakHeightFrac = 0.22;
// Jitter-grid cell size in crater diameters, and a crater's farthest reach in cells.
constexpr double kCellDiams = 2.5;
constexpr double kReachCells = kEjectaQ / (2.0 * kCellDiams);
// Octave fade bounds in vertices across a crater, gone below the low one and full above the high.
constexpr double kFadeLoVerts = 2.0, kFadeHiVerts = 4.0;
// Scatter pass row bands, this many per thread for dynamic balance, never thinner than the floor.
constexpr int kBandsPerThread = 2, kMinBandRows = 16;

// Crater profile at radiusFraction = distance / rim radius. Floor term in [-1,0], rim term in [0,1]. Branch-free so the footprint loop vectorizes.
struct CraterProfile {
    double bowl, rim;
};
inline CraterProfile craterProfile(double radiusFraction, double flatness, double degradation, double rimWidthFraction) {
    const double parabolicBowl = -(1.0 - radiusFraction * radiusFraction);
    constexpr double wallStart = 0.55, invWall = 1.0 / (1.0 - wallStart);
    const double wallFraction = (radiusFraction - wallStart) * invWall;
    const double flatBottomBowl = (radiusFraction < wallStart) ? -1.0 : -(1.0 - wallFraction * wallFraction);
    double freshBowl = parabolicBowl * (1.0 - flatness) + flatBottomBowl * flatness;
    freshBowl = freshBowl * std::abs(freshBowl);   // ease toward the crest
    // Degraded craters lose the sharp wall/rim junction.
    const double degradedBowl = -bump(radiusFraction);
    const double bowl = (radiusFraction < 1.0) ? freshBowl * (1.0 - degradation) + degradedBowl * degradation : 0.0;
    const double rimBump = bump((radiusFraction - 1.0) / (kRimWidth * rimWidthFraction));   // crest at the rim radius
    // Ejecta apron, flat inside the rim and gone at kEjectaQ.
    const double ejectaDistanceFraction = kEjectaQ - radiusFraction;
    const double ejectaEnvelope = (radiusFraction > 1.0) ? ejectaDistanceFraction * ejectaDistanceFraction * (1.0 - 0.6 * degradation * (1.0 - ejectaDistanceFraction)) : 1.0;
    return {bowl, rimBump * ejectaEnvelope};
}

constexpr double kDiamKm[] = {80.0, 40.0, 20.0, 10.0, 5.0, 2.5, 1.25, 0.6, 0.3, 0.15, 0.08, 0.04, 0.02, 0.01, 0.005, 0.0025, 0.00125};
// Occupancy of one cell per class, one meaning a crater in every cell.
constexpr double kDensity[] = {0.275, 0.25, 0.24, 0.24, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
                               0.45, 0.45, 0.6, 0.7, 0.7};
constexpr int kOct = sizeof(kDiamKm) / sizeof(kDiamKm[0]);

static_assert(sizeof(kDensity) / sizeof(kDensity[0]) == kOct,
              "kDensity must have one entry per kDiamKm octave");

struct OctaveInfo {
    double diameterDeg, cellSizeDeg, invCellSizeDeg;
};
inline const OctaveInfo& octave(int octaveIndex) {
    static const auto octaveTable = []() {
        std::array<OctaveInfo, kOct> t{};
        for (int i = 0; i < kOct; ++i) {
            t[i].diameterDeg = kDiamKm[i] / kKmPerDeg;
            t[i].cellSizeDeg = t[i].diameterDeg * kCellDiams;
            t[i].invCellSizeDeg = 1.0 / t[i].cellSizeDeg;
        }
        return t;
    }();
    return octaveTable[octaveIndex];
}
// Octave fade, continuous so neighboring LODs agree.
inline double octaveFade(int octaveIndex, double invSpacingDeg) {
    return qtfield::resolutionFade(octave(octaveIndex).diameterDeg * invSpacingDeg, kFadeLoVerts, kFadeHiVerts);
}

struct CraterInstance {
    double centerLonDeg, centerLatDeg;   // center (deg)
    double nominalRadiusKm;              // nominal radius (km)
    double age, preservationFactor;      // fresh to subdued, and the depth/rim multiplier
    double depthKm, rimKm, flatness;
    bool isComplex;
    double outlineAmplitude, streakAmplitude, invPeakRadiusFraction;
    double cosineHarmonics[6], sineHarmonics[6];   // cos/sin coefficients: outline harmonics 2/3/5, then streak harmonics 5/7/8
    double baseRimWidthFraction;
    double reachKm;   // relief is exactly zero beyond this radius

    // Returns nullopt if the cell has no crater.
    [[nodiscard]] static std::optional<CraterInstance> make(int gx, int gy, int octaveIndex) {
        CraterInstance crater;
        if (hash01(gx, gy, octaveIndex * 13 + 1) > kDensity[octaveIndex]) {
            return std::nullopt;
        }
        const double cellSizeDeg = octave(octaveIndex).cellSizeDeg, classDiameterKm = kDiamKm[octaveIndex];
        // Center anywhere in the cell, or the lattice shows as dimple rows.
        crater.centerLonDeg = (gx + hash01(gx, gy, octaveIndex * 13 + 2)) * cellSizeDeg;
        crater.centerLatDeg = (gy + hash01(gx, gy, octaveIndex * 13 + 3)) * cellSizeDeg;
        crater.nominalRadiusKm = 0.5 * classDiameterKm * (0.55 + 0.45 * hash01(gx, gy, octaveIndex * 13 + 4));
        const double ageSample = hash01(gx, gy, octaveIndex * 13 + 5);
        crater.age = 1.0 - (1.0 - ageSample) * (1.0 - ageSample);   // skewed old
        crater.preservationFactor = 1.0 - 0.85 * crater.age;
        // Depth and rim height from Pike (1977).
        const double diameterKm = 2.0 * crater.nominalRadiusKm;
        crater.isComplex = diameterKm >= 15.0;
        crater.depthKm = (crater.isComplex ? 1.044 * std::pow(diameterKm, 0.301) : 0.196 * diameterKm) * crater.preservationFactor;
        crater.rimKm = (crater.isComplex ? 0.236 * std::pow(diameterKm, 0.399) : 0.036 * diameterKm) * crater.preservationFactor;
        crater.flatness = crater.isComplex ? std::min(1.0, (diameterKm - 15.0) / 30.0) : 0.0;
        crater.outlineAmplitude = 0.04 + 0.08 * crater.age;                              // lumpy outline, old craters more so
        crater.streakAmplitude = (crater.age < 0.8) ? 0.45 * (1.0 - crater.age) : 0.0;   // ejecta rays, fresh only
        const double peakRadiusFraction = kPeakDiamFrac * (1.0 + 0.3 * (hash01(gx, gy, octaveIndex * 13 + 12) - 0.5));
        crater.invPeakRadiusFraction = 1.0 / (peakRadiusFraction * 1.3);
        const double harmonicWeights[3] = {0.50, 0.32, 0.18};
        double outlineMax = 1.0;
        for (int k = 0; k < 3; ++k) {
            crater.cosineHarmonics[k] = 1.4 * harmonicWeights[k] * (hash01(gx, gy, octaveIndex * 13 + 6 + 2 * k) - 0.5);
            crater.sineHarmonics[k] = 1.4 * harmonicWeights[k] * (hash01(gx, gy, octaveIndex * 13 + 7 + 2 * k) - 0.5);
            crater.cosineHarmonics[3 + k] = 1.4 * harmonicWeights[k] * (hash01(gx, gy, octaveIndex * 13 + 14 + 2 * k) - 0.5);
            crater.sineHarmonics[3 + k] = 1.4 * harmonicWeights[k] * (hash01(gx, gy, octaveIndex * 13 + 15 + 2 * k) - 0.5);
            // a cos + b sin <= sqrt(a^2 + b^2), the outline's largest excursion.
            outlineMax += crater.outlineAmplitude * std::sqrt(crater.cosineHarmonics[k] * crater.cosineHarmonics[k] + crater.sineHarmonics[k] * crater.sineHarmonics[k]);
        }
        // Rim half-width, widened per sample so it never aliases.
        crater.baseRimWidthFraction = 0.08 + 0.17 * crater.age;
        crater.reachKm = std::min(kEjectaQ * crater.nominalRadiusKm, kReachQ * crater.nominalRadiusKm * outlineMax);
        return crater;
    }

    [[nodiscard]] double reliefKmAtOffset(double offsetEastKm, double offsetNorthKm, double spacingKm) const {
        const double distanceSquaredKm2 = offsetEastKm * offsetEastKm + offsetNorthKm * offsetNorthKm;
        const double distanceKm = std::sqrt(distanceSquaredKm2);
        const double invDistanceKm = 1.0 / std::max(distanceKm, 1e-9);
        const double cosAzimuth1 = offsetEastKm * invDistanceKm, sinAzimuth1 = offsetNorthKm * invDistanceKm;   // (cos, sin) azimuth
        // Multiple angles by products of lower ones, no trig per sample.
        const double cosAzimuth2 = cosAzimuth1 * cosAzimuth1 - sinAzimuth1 * sinAzimuth1, sinAzimuth2 = 2.0 * cosAzimuth1 * sinAzimuth1;
        const double cosAzimuth3 = cosAzimuth2 * cosAzimuth1 - sinAzimuth2 * sinAzimuth1, sinAzimuth3 = sinAzimuth2 * cosAzimuth1 + cosAzimuth2 * sinAzimuth1;
        const double cosAzimuth5 = cosAzimuth2 * cosAzimuth3 - sinAzimuth2 * sinAzimuth3, sinAzimuth5 = sinAzimuth2 * cosAzimuth3 + cosAzimuth2 * sinAzimuth3;
        const double outline = 1.0 + outlineAmplitude * (cosineHarmonics[0] * cosAzimuth2 + sineHarmonics[0] * sinAzimuth2 + cosineHarmonics[1] * cosAzimuth3 + sineHarmonics[1] * sinAzimuth3 + cosineHarmonics[2] * cosAzimuth5 + sineHarmonics[2] * sinAzimuth5);
        const double invOutlineRadiusKm = 1.0 / (nominalRadiusKm * outline);
        const double radiusFraction = distanceKm * invOutlineRadiusKm;
        const double rimWidthFraction = std::min(kRimSigmaMax, std::max(baseRimWidthFraction, 0.8 * spacingKm * invOutlineRadiusKm));
        const auto [bowl, rim] = craterProfile(radiusFraction, flatness, age, rimWidthFraction);
        // Ejecta rays, past the rim only. Unconditional so the footprint loop vectorizes.
        const double cosAzimuth7 = cosAzimuth2 * cosAzimuth5 - sinAzimuth2 * sinAzimuth5, sinAzimuth7 = sinAzimuth2 * cosAzimuth5 + cosAzimuth2 * sinAzimuth5;
        const double cosAzimuth8 = cosAzimuth3 * cosAzimuth5 - sinAzimuth3 * sinAzimuth5, sinAzimuth8 = sinAzimuth3 * cosAzimuth5 + cosAzimuth3 * sinAzimuth5;
        const double streak = cosineHarmonics[3] * cosAzimuth5 + sineHarmonics[3] * sinAzimuth5 + cosineHarmonics[4] * cosAzimuth7 + sineHarmonics[4] * sinAzimuth7 + cosineHarmonics[5] * cosAzimuth8 + sineHarmonics[5] * sinAzimuth8;
        const double ejectaMod = 1.0 + streakAmplitude * streak * std::min(1.0, std::max(0.0, (radiusFraction - 1.0) * 4.0));
        double reliefKm = depthKm * bowl + rimKm * rim * ejectaMod;
        // Central peak, zero for simple craters and past the peak radius.
        reliefKm += depthKm * kPeakHeightFrac * preservationFactor * flatness * bump(radiusFraction * invPeakRadiusFraction);
        return reliefKm;
    }
};

}   // namespace

// Reach is at most kReachCells of a cell, so only adjacent cells can contribute.
double CraterField::heightAt(double lonDeg, double latDeg, double sampleSpacingDeg) {
    const double invSampleSpacingDeg = 1.0 / std::max(sampleSpacingDeg, 1e-9);
    const double cosLat = qtplanet::cosLatClamped(latDeg);
    const double spacingKm = sampleSpacingDeg * kKmPerDeg * cosLat;
    double reliefM = 0.0;
    for (int octaveIndex = 0; octaveIndex < kOct; ++octaveIndex) {
        const double fade = octaveFade(octaveIndex, invSampleSpacingDeg);
        if (fade <= 0.0) {
            break;
        }
        const double cx = lonDeg * octave(octaveIndex).invCellSizeDeg, cy = latDeg * octave(octaveIndex).invCellSizeDeg;
        const int ix = fastFloor(cx), iy = fastFloor(cy);
        const double fx = cx - ix, fy = cy - iy;
        // Walk the neighboring cells, but skip ones out of reach.
        const int ox0 = fx < kReachCells ? -1 : 0, ox1 = fx > 1.0 - kReachCells ? 1 : 0;
        const int oy0 = fy < kReachCells ? -1 : 0, oy1 = fy > 1.0 - kReachCells ? 1 : 0;
        for (int oy = oy0; oy <= oy1; ++oy) {
            for (int ox = ox0; ox <= ox1; ++ox) {
                const auto candidate = CraterInstance::make(ix + ox, iy + oy, octaveIndex);
                if (!candidate) {
                    continue;
                }
                const CraterInstance& crater = *candidate;
                const double dLon = qtplanet::wrapLongitude(lonDeg - crater.centerLonDeg);
                const double offsetEastKm = dLon * cosLat * kKmPerDeg, offsetNorthKm = (latDeg - crater.centerLatDeg) * kKmPerDeg;
                if (offsetEastKm * offsetEastKm + offsetNorthKm * offsetNorthKm < crater.reachKm * crater.reachKm) {
                    reliefM += crater.reliefKmAtOffset(offsetEastKm, offsetNorthKm, spacingKm) * 1000.0 * fade;
                }
            }
        }
    }
    reliefM += RockField::beddingAt(lonDeg, latDeg, sampleSpacingDeg);
    reliefM += SurfaceRoughness::heightAt(lonDeg, latDeg, sampleSpacingDeg);
    return reliefM;
}

namespace {

// Row bands for the three scatter passes, sized by kBandsPerThread and kMinBandRows.
struct RowBands {
    int sampleCount = 0, rows = 0;
    int count() const { return rows > 0 ? (sampleCount + rows - 1) / rows : 0; }
    int begin(int band) const { return band * rows; }
    int end(int band) const { return std::min(sampleCount, (band + 1) * rows); }
};

RowBands rowBands(int sampleCount) {
    if (sampleCount <= 0) {
        return {0, 0};
    }
    if (!qtplanet::parallelGrid(static_cast<size_t>(sampleCount) * sampleCount)) {
        return {sampleCount, sampleCount};
    }
    const int bands = kBandsPerThread * qtplanet::parallelThreads();
    return {sampleCount, std::max(kMinBandRows, (sampleCount + bands - 1) / bands)};
}

// One octave's craters over the whole grid, a vector per cell row in gx order, the order heightAt() sums in.
struct OctaveCraters {
    double scale, cellSizeDeg, invCellSizeDeg;
    int gy0, gy1;
    std::vector<std::vector<CraterInstance>> rows;
};

// Instantiates every crater that can reach the grid once, cell rows in parallel. The bands then walk their own rows.
std::vector<OctaveCraters> instantiate(double originLonDeg, double originLatDeg, double lon1, double lat1, double invSampleSpacingDeg) {
    std::vector<OctaveCraters> out;
    std::vector<std::pair<int, int>> jobs;   // (octave, row)
    for (int octaveIndex = 0; octaveIndex < kOct; ++octaveIndex) {
        const double fade = octaveFade(octaveIndex, invSampleSpacingDeg);
        if (fade <= 0.0) {
            break;
        }
        OctaveCraters octaveInfo;
        octaveInfo.scale = 1000.0 * fade;
        octaveInfo.cellSizeDeg = octave(octaveIndex).cellSizeDeg;
        octaveInfo.invCellSizeDeg = octave(octaveIndex).invCellSizeDeg;
        octaveInfo.gy0 = fastFloor((originLatDeg - kReachCells * octaveInfo.cellSizeDeg) * octaveInfo.invCellSizeDeg);
        octaveInfo.gy1 = fastFloor((lat1 + kReachCells * octaveInfo.cellSizeDeg) * octaveInfo.invCellSizeDeg);
        octaveInfo.rows.resize(static_cast<size_t>(octaveInfo.gy1 - octaveInfo.gy0 + 1));
        for (int gy = octaveInfo.gy0; gy <= octaveInfo.gy1; ++gy) {
            jobs.emplace_back(octaveIndex, gy);
        }
        out.push_back(std::move(octaveInfo));
    }
    const int jobCount = static_cast<int>(jobs.size());
#pragma omp parallel for num_threads(qtplanet::parallelThreads()) schedule(dynamic, 4) if (jobCount >= 64)
    for (int t = 0; t < jobCount; ++t) {
        const int octaveIndex = jobs[t].first, gy = jobs[t].second;
        OctaveCraters& octaveInfo = out[octaveIndex];
        const int gx0 = fastFloor((originLonDeg - kReachCells * octaveInfo.cellSizeDeg) * octaveInfo.invCellSizeDeg), gx1 = fastFloor((lon1 + kReachCells * octaveInfo.cellSizeDeg) * octaveInfo.invCellSizeDeg);
        std::vector<CraterInstance>& row = octaveInfo.rows[gy - octaveInfo.gy0];
        for (int gx = gx0; gx <= gx1; ++gx) {
            if (auto crater = CraterInstance::make(gx, gy, octaveIndex)) {
                row.push_back(*crater);
            }
        }
    }
    return out;
}

// Crater relief over rows [rowBegin, rowEnd) of the grid, craters in the whole-grid walk order, so sums match serial bit for bit.
void addCraterRows(const qtfield::GridSpec& grid, qtfield::RowRange rows,
                   const std::vector<OctaveCraters>& craterOctaves, const double* kmPerDegLon,
                   const double* spacingKm, std::vector<double>& heightsM) {
    const auto [originLonDeg, originLatDeg, stepLonDeg, stepLatDeg, sampleCount, sampleSpacingDeg] = grid;
    const auto [rowBegin, rowEnd] = rows;
    const double latB0 = originLatDeg + rowBegin * stepLatDeg, latB1 = originLatDeg + (rowEnd - 1) * stepLatDeg;
    const double invStepLon = 1.0 / stepLonDeg, invStepLat = 1.0 / stepLatDeg;
    for (const OctaveCraters& octaveInfo : craterOctaves) {
        const double scale = octaveInfo.scale;
        // Cell rows whose crater can touch the band.
        const int gyB0 = std::max(octaveInfo.gy0, fastFloor((latB0 - kReachCells * octaveInfo.cellSizeDeg) * octaveInfo.invCellSizeDeg));
        const int gyB1 = std::min(octaveInfo.gy1, fastFloor((latB1 + kReachCells * octaveInfo.cellSizeDeg) * octaveInfo.invCellSizeDeg));
        for (int gy = gyB0; gy <= gyB1; ++gy) {
            for (const CraterInstance& cref : octaveInfo.rows[gy - octaveInfo.gy0]) {
                const CraterInstance crater = cref;   // a local, since through the vector it could alias the rows and the loop stops vectorizing
                // Walk the disk where the dH is nonzero, row by row.
                const double reach2 = crater.reachKm * crater.reachKm;
                const double reachLat = crater.reachKm / kKmPerDeg;
                const int j0 = std::max(rowBegin, static_cast<int>(std::ceil((crater.centerLatDeg - reachLat - originLatDeg) * invStepLat)));
                const int j1 = std::min(rowEnd - 1, static_cast<int>(std::floor((crater.centerLatDeg + reachLat - originLatDeg) * invStepLat)));
                for (int j = j0; j <= j1; ++j) {
                    const double offsetNorthKm = (originLatDeg + j * stepLatDeg - crater.centerLatDeg) * kKmPerDeg;
                    const double half2 = reach2 - offsetNorthKm * offsetNorthKm;
                    if (half2 <= 0.0) {
                        continue;
                    }
                    const double kmPerDeg = kmPerDegLon[j], sampleSpacingKm = spacingKm[j];
                    const double halfLon = std::sqrt(half2) / kmPerDeg;
                    const int i0 = std::max(0, static_cast<int>(std::ceil((crater.centerLonDeg - halfLon - originLonDeg) * invStepLon)));
                    const int i1 = std::min(sampleCount - 1, static_cast<int>(std::floor((crater.centerLonDeg + halfLon - originLonDeg) * invStepLon)));
                    double* row = &heightsM[static_cast<size_t>(j) * sampleCount];
                    for (int i = i0; i <= i1; ++i) {
                        const double offsetEastKm = (originLonDeg + i * stepLonDeg - crater.centerLonDeg) * kmPerDeg;
                        row[i] += crater.reliefKmAtOffset(offsetEastKm, offsetNorthKm, sampleSpacingKm) * scale;
                    }
                }
            }
        }
    }
}

}   // namespace

// Row bands in parallel. Each band adds craters, then plinths, then roughness, in heightAt() order.
void CraterField::addGrid(const qtfield::GridSpec& grid, std::vector<double>& heightsM) {
    const auto [originLonDeg, originLatDeg, stepLonDeg, stepLatDeg, sampleCount, sampleSpacingDeg] = grid;
    std::vector<double> kmPerDegLon(sampleCount), spacingKm(sampleCount);
    for (int j = 0; j < sampleCount; ++j) {
        const double cosLat = qtplanet::cosLatClamped(originLatDeg + j * stepLatDeg);
        kmPerDegLon[j] = cosLat * kKmPerDeg;
        spacingKm[j] = sampleSpacingDeg * kKmPerDeg * cosLat;
    }
    const auto craterOctaves = instantiate(originLonDeg, originLatDeg, originLonDeg + (sampleCount - 1) * stepLonDeg, originLatDeg + (sampleCount - 1) * stepLatDeg, 1.0 / std::max(sampleSpacingDeg, 1e-9));
    const RowBands bands = rowBands(sampleCount);
    const int bandCount = bands.count();
#pragma omp parallel for num_threads(qtplanet::parallelThreads()) schedule(dynamic) if (bandCount > 1)
    for (int b = 0; b < bandCount; ++b) {
        const int rowBegin = bands.begin(b), rowEnd = bands.end(b);
        addCraterRows(grid, {rowBegin, rowEnd}, craterOctaves, kmPerDegLon.data(), spacingKm.data(), heightsM);
        RockField::addBeddingGridRows(grid, {rowBegin, rowEnd}, heightsM);
        SurfaceRoughness::addGridRows(grid, {rowBegin, rowEnd}, heightsM);
    }
}
