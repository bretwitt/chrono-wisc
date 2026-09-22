#include "chrono_planet/dem/MultiGeoTIFFManager.h"

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "chrono_planet/core/MathUtil.h"
#include "chrono_planet/core/Parallel.h"
#include "chrono_planet/core/Planet.h"

namespace {

// Elevations below this are treated as nodata whatever the band declares.
constexpr double kInvalidBelowM = -30000.0;
// Overlay blend ramp on the ground, and its bounds in pixels: the B-spline support at least, one byte at most.
constexpr double kFeatherM = 500.0;
constexpr double kMinFeatherPx = 4.0, kMaxFeatherPx = 255.0;
// Chamfer distance transform steps in thirds of a pixel, and the value of a pixel it has not reached.
constexpr int kChamferAxial = 3, kChamferDiagonal = 4, kChamferUnreached = 100000;
// The B-spline footprint is 4x4, so a sample this far inside valid data is fully trusted.
constexpr double kFullSupportPx = 2.0;
// Geotransforms with a smaller determinant are degenerate.
constexpr double kMinGeoTransformDet = 1e-12;
// Raster axis skew below this counts as axis-aligned.
constexpr double kMaxAxisSkew = 1e-12;

// Cubic B-spline weights of the four taps around a sample at fractional offset t in [0, 1).
std::array<double, 4> bsplineWeights(double t) {
    const double t2 = t * t, t3 = t2 * t;
    return {(1.0 - 3.0 * t + 3.0 * t2 - t3) / 6.0, (4.0 - 6.0 * t2 + 3.0 * t3) / 6.0,
            (1.0 + 3.0 * t + 3.0 * t2 - 3.0 * t3) / 6.0, t3 / 6.0};
}

// Raster (col, row) of a dataset-space point.
struct GeoTransformInverse {
    std::array<double, 6> gt;
    double det;
    explicit GeoTransformInverse(const std::array<double, 6>& g) : gt(g), det(g[1] * g[5] - g[2] * g[4]) {}
    bool valid() const { return std::abs(det) >= kMinGeoTransformDet; }
    double col(double x, double y) const { return ((x - gt[0]) * gt[5] - (y - gt[3]) * gt[2]) / det; }
    double row(double x, double y) const { return (-(x - gt[0]) * gt[4] + (y - gt[3]) * gt[1]) / det; }
};

// Reports once a source whose raster axes are not aligned with lon/lat. It is treated as not covering.
void warnNonseparableRaster() {
    static std::once_flag warned;
    std::call_once(warned, [] {
        std::cerr << "[MGTM] raster is not separable in lon/lat; reproject it to Plate Carree\n";
    });
}

}   // namespace

MultiGeoTIFFManager::MultiGeoTIFFManager() {
    moonLonLat_.SetGeogCS("GCS_Moon_2000", "D_Moon_2000", "Moon_2000_IAU_IAG", qtplanet::kRadiusM, 0.0);
    moonLonLat_.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
}

bool MultiGeoTIFFManager::addSource(const std::string& file, int minZoom, int maxZoom) {
    auto ldr = std::make_shared<GeoTIFFLoader>();
    try {
        ldr->load(file);
    } catch (const std::exception& e) {
        std::cerr << "[MGTM] load failed: " << e.what() << '\n';
        return false;
    }
    const OGRSpatialReference& dsSRS = ldr->srs();
    const bool dsIsGeo = dsSRS.IsGeographic();

    TransformPtr toDataset, toLonLat;
    if (!dsIsGeo) {
        toDataset.reset(OGRCreateCoordinateTransformation(&moonLonLat_, &dsSRS));
        toLonLat.reset(OGRCreateCoordinateTransformation(&dsSRS, &moonLonLat_));
        if (!toDataset || !toLonLat) {
            std::cerr << "[MGTM] can't build transform into " << file << '\n';
            return false;
        }
    }

    // Lon/lat extent from the raster corners.
    const auto& gt = ldr->geoTransform();
    const double w = ldr->width() - 1, h = ldr->height() - 1;
    double lon[4] = {gt[0], gt[0] + w * gt[1], gt[0] + h * gt[2], gt[0] + w * gt[1] + h * gt[2]};
    double lat[4] = {gt[3], gt[3] + w * gt[4], gt[3] + h * gt[5], gt[3] + w * gt[4] + h * gt[5]};
    if (toLonLat) {
        toLonLat->Transform(4, lon, lat);
    }

    Source src{ldr, minZoom, maxZoom,
               *std::min_element(lon, lon + 4), *std::max_element(lon, lon + 4),
               *std::min_element(lat, lat + 4), *std::max_element(lat, lat + 4),
               dsIsGeo, std::move(toDataset)};
    // The base source needs no mask. Overlays fade to it at edges and holes.
    if (minZoom > 0) {
        const double pixelMetres = std::abs(gt[1]) * (dsIsGeo ? qtplanet::kMetresPerDegree : 1.0);
        src.mask = validityMaskOf(*ldr, pixelMetres);
    }
    sources_.push_back(std::move(src));

    std::cout << "[MGTM] added source: " << file
              << (ldr->hasNodata() ? "  (nodata " + std::to_string(ldr->nodata()) + ")" : "") << '\n';
    return true;
}

// Chamfer distance transform to the nearest invalid pixel or raster edge, in whole pixels.
MultiGeoTIFFManager::ValidityMask MultiGeoTIFFManager::validityMaskOf(const GeoTIFFLoader& raster, double pixelMetres) {
    const int w = raster.width(), h = raster.height();
    const auto& d = raster.elevation();
    const bool hasNodata = raster.hasNodata();
    const double nodata = raster.nodata();
    ValidityMask mask;
    mask.featherPx = static_cast<float>(std::clamp(kFeatherM / std::max(pixelMetres, 1e-3), kMinFeatherPx, kMaxFeatherPx));

    std::vector<int> dist(static_cast<size_t>(w) * h);   // thirds of a pixel
    auto at = [&](int x, int y) -> int& { return dist[static_cast<size_t>(y) * w + x]; };
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            const double v = d[static_cast<size_t>(y) * w + x];
            const bool invalid = !std::isfinite(v) || (hasNodata && v == nodata) || v < kInvalidBelowM;
            const bool edge = x == 0 || y == 0 || x == w - 1 || y == h - 1;
            at(x, y) = invalid ? 0 : edge ? kChamferAxial
                                          : kChamferUnreached;
        }
    }
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            int& v = at(x, y);
            if (x > 0) {
                v = std::min(v, at(x - 1, y) + kChamferAxial);
            }
            if (y > 0) {
                v = std::min(v, at(x, y - 1) + kChamferAxial);
            }
            if (x > 0 && y > 0) {
                v = std::min(v, at(x - 1, y - 1) + kChamferDiagonal);
            }
            if (x + 1 < w && y > 0) {
                v = std::min(v, at(x + 1, y - 1) + kChamferDiagonal);
            }
        }
    }
    for (int y = h - 1; y >= 0; --y) {
        for (int x = w - 1; x >= 0; --x) {
            int& v = at(x, y);
            if (x + 1 < w) {
                v = std::min(v, at(x + 1, y) + kChamferAxial);
            }
            if (y + 1 < h) {
                v = std::min(v, at(x, y + 1) + kChamferAxial);
            }
            if (x + 1 < w && y + 1 < h) {
                v = std::min(v, at(x + 1, y + 1) + kChamferDiagonal);
            }
            if (x > 0 && y + 1 < h) {
                v = std::min(v, at(x - 1, y + 1) + kChamferDiagonal);
            }
        }
    }
    mask.distancePx.resize(dist.size());
    size_t invalidCount = 0;
    for (size_t i = 0; i < dist.size(); ++i) {
        mask.distancePx[i] = static_cast<unsigned char>(std::min(255, dist[i] / kChamferAxial));
        if (dist[i] == 0) {
            ++invalidCount;
        }
    }
    std::cout << "[MGTM]   validity mask: " << w << "x" << h << ", " << (100.0 * invalidCount / dist.size())
              << "% nodata, feather " << mask.featherPx << " px\n";
    return mask;
}

double MultiGeoTIFFManager::weightAt(const Source& s, double colF, double rowF) {
    if (s.mask.complete()) {
        return 1.0;
    }
    const int w = s.loader->width(), h = s.loader->height();
    if (!(colF >= 0.0) || !(rowF >= 0.0) || colF >= w - 1 || rowF >= h - 1) {
        return 0.0;
    }
    const int c0 = static_cast<int>(colF), r0 = static_cast<int>(rowF);
    const double tx = colF - c0, ty = rowF - r0;
    auto v = [&](int c, int r) { return static_cast<double>(s.mask.distancePx[static_cast<size_t>(r) * w + c]); };
    const double dpx = (1 - tx) * (1 - ty) * v(c0, r0) + tx * (1 - ty) * v(c0 + 1, r0) +
                       (1 - tx) * ty * v(c0, r0 + 1) + tx * ty * v(c0 + 1, r0 + 1);
    // Fully trusted once the whole B-spline footprint is on valid data, then ramp.
    const double t = std::clamp((dpx - kFullSupportPx) / std::max(static_cast<double>(s.mask.featherPx) - kFullSupportPx, 1.0), 0.0, 1.0);
    return qtplanet::smoothstep01(t);
}

MultiGeoTIFFManager::SelectedSources MultiGeoTIFFManager::selectSources(double minLon, double maxLon, double minLat, double maxLat, int zoom) const {
    SelectedSources selected;
    auto& [fine, base] = selected;
    int bestFineSpan = std::numeric_limits<int>::max(), bestBaseMax = std::numeric_limits<int>::min();
    for (const auto& s : sources_) {
        bool intersects = false, covers = false;
        if (!(maxLat < s.minLat || minLat > s.maxLat)) {
            for (double shift : {0.0, 360.0, -360.0}) {
                const double a = minLon + shift, b = maxLon + shift;
                if (!(b < s.minLon || a > s.maxLon)) {
                    intersects = true;
                }
                if (a >= s.minLon && b <= s.maxLon && minLat >= s.minLat && maxLat <= s.maxLat) {
                    covers = true;
                }
            }
        }
        if (zoom >= s.minZoom && zoom <= s.maxZoom) {
            if (intersects && s.maxZoom - s.minZoom < bestFineSpan) {
                bestFineSpan = s.maxZoom - s.minZoom;
                fine = &s;
            }
        } else if (s.maxZoom < zoom && covers && s.maxZoom > bestBaseMax) {
            bestBaseMax = s.maxZoom;
            base = &s;
        }
    }
    // A complete fine source has nothing to fade to, so treat it as base and sample once.
    if (fine && fine->mask.complete()) {
        base = fine;
        fine = nullptr;
    }
    return selected;
}

// Separable B-spline reconstruction. Projected sources are checked for separability on the grid's corners.
std::optional<MultiGeoTIFFManager::WeightedGrid> MultiGeoTIFFManager::gridReconstruct(const Source& src, double lon0, double lat0, double stepLon, double stepLat,
                                                                                      int nLon, int nLat) const {
    WeightedGrid grid;
    auto& out = grid.elevations;
    auto& wt = grid.weights;
    const GeoTransformInverse inv(src.loader->geoTransform());
    if (!inv.valid()) {
        return std::nullopt;
    }

    double shift = 0.0;
    if (lon0 < src.minLon) {
        shift = 360.0;
    } else if (lon0 > src.maxLon) {
        shift = -360.0;
    }
    std::vector<double> cx(nLon), cy(nLat);
    for (int i = 0; i < nLon; ++i) {
        cx[i] = lon0 + shift + i * stepLon;
    }
    for (int j = 0; j < nLat; ++j) {
        cy[j] = lat0 + j * stepLat;
    }

    // Raster x of every grid column, raster y of every grid row.
    std::vector<double> colF(nLon), rowF(nLat);
    if (src.isGeographic) {
        const auto& gt = inv.gt;
        if (std::abs(gt[2]) > kMaxAxisSkew || std::abs(gt[4]) > kMaxAxisSkew) {
            warnNonseparableRaster();
            return std::nullopt;   // rotated raster
        }
        for (int i = 0; i < nLon; ++i) {
            colF[i] = inv.col(cx[i], cy[0]);
        }
        for (int j = 0; j < nLat; ++j) {
            rowF[j] = inv.row(cx[0], cy[j]);
        }
    } else {
        // Transform the first row and column, then check the far corners agree.
        std::vector<double> xs(cx), ys(nLon, cy[0]);
        std::vector<double> xs2(nLat, cx[0]), ys2(cy);
        double cxs[2] = {cx[nLon - 1], cx[0]}, cys[2] = {cy[nLat - 1], cy[nLat - 1]};
        {
            std::lock_guard<std::mutex> lock(transformMutex_);
            if (!src.toDataset->Transform(nLon, xs.data(), ys.data())) {
                return std::nullopt;
            }
            if (!src.toDataset->Transform(nLat, xs2.data(), ys2.data())) {
                return std::nullopt;
            }
            if (!src.toDataset->Transform(2, cxs, cys)) {
                return std::nullopt;
            }
        }
        const double tol = 1e-6 * (std::abs(xs[nLon - 1] - xs[0]) + std::abs(ys2[nLat - 1] - ys2[0]) + 1.0);
        const bool separable = std::abs(cxs[0] - xs[nLon - 1]) < tol && std::abs(cys[0] - ys2[nLat - 1]) < tol &&
                               std::abs(cxs[1] - xs[0]) < tol && std::abs(cys[1] - ys2[nLat - 1]) < tol;
        if (!separable) {
            warnNonseparableRaster();
            return std::nullopt;
        }
        for (int i = 0; i < nLon; ++i) {
            colF[i] = inv.col(xs[i], ys[0]);
        }
        for (int j = 0; j < nLat; ++j) {
            rowF[j] = inv.row(xs2[0], ys2[j]);
        }
    }

    const int w = src.loader->width(), h = src.loader->height();
    const auto& d = src.loader->elevation();
    const size_t n = static_cast<size_t>(nLon) * nLat;
    out.assign(n, 0.0);
    wt.assign(n, 0.0);

    // Per-column base index and weights, per-row likewise.
    std::vector<int> c0(nLon), r0(nLat);
    std::vector<double> wx(static_cast<size_t>(nLon) * 4), wy(static_cast<size_t>(nLat) * 4);
    std::vector<unsigned char> colMissing(nLon, 0), rowMissing(nLat, 0);
    int cMin = INT_MAX, cMax = INT_MIN;
    for (int i = 0; i < nLon; ++i) {
        c0[i] = int(std::floor(colF[i]));
        if (c0[i] < 0 || c0[i] + 1 >= w) {
            colMissing[i] = 1;
            continue;
        }
        const std::array<double, 4> weights = bsplineWeights(colF[i] - c0[i]);
        std::copy(weights.begin(), weights.end(), &wx[i * 4]);
        cMin = std::min(cMin, c0[i]);
        cMax = std::max(cMax, c0[i]);
    }
    int rMin = INT_MAX, rMax = INT_MIN;
    for (int j = 0; j < nLat; ++j) {
        r0[j] = int(std::floor(rowF[j]));
        if (r0[j] < 0 || r0[j] + 1 >= h) {
            rowMissing[j] = 1;
            continue;
        }
        const std::array<double, 4> weights = bsplineWeights(rowF[j] - r0[j]);
        std::copy(weights.begin(), weights.end(), &wy[j * 4]);
        rMin = std::min(rMin, r0[j] - 1);
        rMax = std::max(rMax, r0[j] + 2);
    }
    if (cMin > cMax || rMin > rMax) {
        return grid;   // grid entirely outside, all weights 0
    }
    rMin = std::clamp(rMin, 0, h - 1);
    rMax = std::clamp(rMax, 0, h - 1);
    const int nRows = rMax - rMin + 1;

    // Per-column taps, clamped once so the row pass is straight-line.
    std::vector<int> tap(static_cast<size_t>(nLon) * 4);
    for (int i = 0; i < nLon; ++i) {
        for (int k = 0; k < 4; ++k) {
            tap[i * 4 + k] = colMissing[i] ? 0 : std::clamp(c0[i] - 1 + k, 0, w - 1);
        }
    }
    // Row pass, raster rows rMin..rMax filtered along x at every grid column. Both passes are row-independent.
    std::vector<double> tmp(static_cast<size_t>(nRows) * nLon);
    const bool parRows = qtplanet::parallelGrid(tmp.size());
#pragma omp parallel for num_threads(qtplanet::parallelThreads()) schedule(static) if (parRows)
    for (int rr = 0; rr < nRows; ++rr) {
        const double* row = &d[static_cast<size_t>(rMin + rr) * w];
        double* trow = &tmp[static_cast<size_t>(rr) * nLon];
        for (int i = 0; i < nLon; ++i) {
            const double* wtp = &wx[i * 4];
            const int* t = &tap[i * 4];
            double acc = 0.0;
            for (int k = 0; k < 4; ++k) {
                acc += wtp[k] * row[t[k]];
            }
            trow[i] = colMissing[i] ? 0.0 : acc;
        }
    }
    // Column pass, with the validity weight per sample.
    const bool masked = !src.mask.complete();
    const bool parCols = qtplanet::parallelGrid(n);
#pragma omp parallel for num_threads(qtplanet::parallelThreads()) schedule(static) if (parCols)
    for (int j = 0; j < nLat; ++j) {
        if (rowMissing[j]) {
            continue;
        }
        const double* wtp = &wy[j * 4];
        const double* rows[4];
        for (int m = 0; m < 4; ++m) {
            rows[m] = &tmp[static_cast<size_t>(std::clamp(r0[j] - 1 + m, 0, h - 1) - rMin) * nLon];
        }
        double* orow = &out[static_cast<size_t>(j) * nLon];
        double* wrow = &wt[static_cast<size_t>(j) * nLon];
        // Weights first, then the straight-line combine keyed on them.
        for (int i = 0; i < nLon; ++i) {
            wrow[i] = 0.0;
            if (!colMissing[i]) {
                wrow[i] = masked ? weightAt(src, colF[i], rowF[j]) : 1.0;
            }
        }
        for (int i = 0; i < nLon; ++i) {
            const double v = wtp[0] * rows[0][i] + wtp[1] * rows[1][i] + wtp[2] * rows[2][i] + wtp[3] * rows[3][i];
            const bool on = wrow[i] > 0.0;
            orow[i] = on ? v : 0.0;
            wrow[i] = on ? wrow[i] : 0.0;
        }
    }
    return grid;
}

std::optional<MultiGeoTIFFManager::GridSamples> MultiGeoTIFFManager::sampleGrid(
    double lon0, double lat0, double stepLon, double stepLat, int nLon, int nLat, int zoom) const {
    if (nLon <= 0 || nLat <= 0) {
        throw std::invalid_argument("sampleGrid requires positive dimensions");
    }
    const double lon1 = lon0 + (nLon - 1) * stepLon, lat1 = lat0 + (nLat - 1) * stepLat;
    const auto [fine, base] = selectSources(std::min(lon0, lon1), std::max(lon0, lon1),
                                            std::min(lat0, lat1), std::max(lat0, lat1), zoom);
    if (!fine && !base) {
        return std::nullopt;
    }
    const size_t n = static_cast<size_t>(nLon) * nLat;
    WeightedGrid blended;
    if (base) {
        auto grid = gridReconstruct(*base, lon0, lat0, stepLon, stepLat, nLon, nLat);
        if (!grid) {
            return std::nullopt;
        }
        blended = std::move(*grid);
    } else {
        blended.elevations.assign(n, 0.0);
        blended.weights.assign(n, 0.0);
    }
    if (fine) {
        const auto grid = gridReconstruct(*fine, lon0, lat0, stepLon, stepLat, nLon, nLat);
        if (!grid && !base) {
            return std::nullopt;
        }
        if (grid) {
            for (size_t i = 0; i < n; ++i) {
                if (grid->weights[i] <= 0.0) {
                    continue;
                }
                const double weight = blended.weights[i] > 0.0 ? grid->weights[i] : 1.0;
                blended.elevations[i] = blended.elevations[i] * (1.0 - weight) + grid->elevations[i] * weight;
                blended.weights[i] = std::max(blended.weights[i], grid->weights[i]);
            }
        }
    }
    // Even if the fine source failed, report holes in the retained base grid.
    GridSamples result;
    result.elevations = std::move(blended.elevations);
    for (size_t i = 0; i < n; ++i) {
        if (blended.weights[i] <= 0.0) {
            result.missing.push_back(i);
        }
    }
    return result;
}

std::optional<double> MultiGeoTIFFManager::sample(double lonDeg, double latDeg, int zoom) const {
    const auto grid = sampleGrid(lonDeg, latDeg, 1e-7, 1e-7, 1, 1, zoom);
    if (!grid || !grid->missing.empty()) {
        return std::nullopt;
    }
    return grid->elevations[0];
}
