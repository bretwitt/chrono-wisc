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
// Test helpers for the Planet module: synthetic GeoTIFFs written with GDAL,
// and the fixtures of the Moon regression test.
//
// =============================================================================

#ifndef PLANET_TEST_UTILS_H
#define PLANET_TEST_UTILS_H

#include <cmath>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include <unistd.h>

namespace planet_test {

constexpr double kPi = 3.14159265358979323846;

/// A directory for this process's test rasters, created on first use.
inline std::string TempDir() {
    static const std::string dir = [] {
        auto p = std::filesystem::temp_directory_path() / ("chrono_planet_utest_" + std::to_string(::getpid()));
        std::filesystem::create_directories(p);
        return p.string();
    }();
    return dir;
}

/// Options of a raster written by WriteRaster.
struct RasterOptions {
    std::optional<double> nodata;
    std::optional<double> band_scale;
    std::optional<double> band_offset;
    int bands = 1;  ///< the data goes into the last band; earlier bands hold zeros
};

/// Write a Float32 GeoTIFF with the given geotransform and coordinate system (null: none).
inline void WriteRaster(const std::string& path,
                        int w,
                        int h,
                        const double gt[6],
                        const OGRSpatialReference* srs,
                        const std::vector<float>& data,
                        const RasterOptions& opt = {}) {
    GDALAllRegister();
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("GTiff");
    GDALDataset* ds = drv->Create(path.c_str(), w, h, opt.bands, GDT_Float32, nullptr);
    if (!ds)
        throw std::runtime_error("cannot create " + path);
    ds->SetGeoTransform(const_cast<double*>(gt));
    if (srs)
        ds->SetSpatialRef(srs);
    std::vector<float> zeros(static_cast<size_t>(w) * h, 0.f);
    for (int b = 1; b <= opt.bands; ++b) {
        GDALRasterBand* band = ds->GetRasterBand(b);
        const bool last = b == opt.bands;
        if (last && opt.nodata)
            band->SetNoDataValue(*opt.nodata);
        if (last && opt.band_scale)
            band->SetScale(*opt.band_scale);
        if (last && opt.band_offset)
            band->SetOffset(*opt.band_offset);
        float* src = last ? const_cast<float*>(data.data()) : zeros.data();
        if (band->RasterIO(GF_Write, 0, 0, w, h, src, w, h, GDT_Float32, 0, 0) != CE_None)
            throw std::runtime_error("cannot write " + path);
    }
    GDALClose(ds);
}

/// A geographic coordinate system on a sphere of the given radius.
inline OGRSpatialReference SphereGeog(double radius) {
    OGRSpatialReference geo;
    geo.SetGeogCS("GCS_Test", "D_Test", "Test_Sphere", radius, 0.0);
    return geo;
}

/// Write a global geographic raster at `res` degrees per pixel, starting at longitude lon_start,
/// holding f(lon, lat) with an optional transform of the stored value.
template <class F>
void WriteGlobal(const std::string& path,
                 double radius,
                 double res,
                 double lon_start,
                 F f,
                 const RasterOptions& opt = {},
                 bool with_srs = true) {
    const int w = static_cast<int>(std::round(360.0 / res)), h = static_cast<int>(std::round(180.0 / res));
    const double gt[6] = {lon_start, res, 0.0, 90.0, 0.0, -res};
    std::vector<float> data(static_cast<size_t>(w) * h);
    for (int r = 0; r < h; ++r)
        for (int c = 0; c < w; ++c)
            data[static_cast<size_t>(r) * w + c] = static_cast<float>(f(lon_start + c * res, 90.0 - r * res));
    const OGRSpatialReference srs = SphereGeog(radius);
    WriteRaster(path, w, h, gt, with_srs ? &srs : nullptr, data, opt);
}

// -----------------------------------------------------------------------------
// Moon regression fixtures. Changing any of these invalidates the reference values.

namespace golden {

constexpr double kMoonRadius = 1737400.0;
constexpr double kSiteLon = 30.75, kSiteLat = 20.19;
constexpr double kOverlayPixel = 50.0;
constexpr int kOverlaySide = 400;

inline double BaseHeight(double lon_deg, double lat_deg) {
    const double lon = lon_deg * kPi / 180, lat = lat_deg * kPi / 180;
    return 1500.0 * std::sin(2 * lon) * std::cos(lat) + 800.0 * std::cos(3 * lat);
}

/// Global geographic base at 1 degree, and a projected 20 km overlay around the site with a nodata hole.
inline void WriteDems(const std::string& base, const std::string& overlay, double radius) {
    const OGRSpatialReference geo = SphereGeog(radius);
    {
        const int w = 360, h = 180;
        const double gt[6] = {-180.0, 1.0, 0.0, 90.0, 0.0, -1.0};
        std::vector<float> data(static_cast<size_t>(w) * h);
        for (int r = 0; r < h; ++r)
            for (int c = 0; c < w; ++c)
                data[static_cast<size_t>(r) * w + c] = static_cast<float>(BaseHeight(-180.0 + c, 90.0 - r));
        WriteRaster(base, w, h, gt, &geo, data);
    }
    {
        OGRSpatialReference proj;
        proj.SetProjCS("Test_Equirectangular");
        proj.SetGeogCS("GCS_Test", "D_Test", "Test_Sphere", radius, 0.0);
        proj.SetEquirectangular2(0.0, 0.0, 0.0, 0.0, 0.0);
        const int n = kOverlaySide;
        const double cx = radius * kSiteLon * kPi / 180, cy = radius * kSiteLat * kPi / 180;
        const double half = 0.5 * n * kOverlayPixel;
        const double gt[6] = {cx - half, kOverlayPixel, 0.0, cy + half, 0.0, -kOverlayPixel};
        std::vector<float> data(static_cast<size_t>(n) * n);
        for (int r = 0; r < n; ++r)
            for (int c = 0; c < n; ++c) {
                const double x = gt[0] + (c + 0.5) * kOverlayPixel, y = gt[3] - (r + 0.5) * kOverlayPixel;
                const double lon = x / radius * 180 / kPi, lat = y / radius * 180 / kPi;
                double v = BaseHeight(lon, lat) + 30.0 * std::sin(x / 700.0) * std::cos(y / 900.0);
                const double dc = c - 0.8 * n, dr = r - 0.25 * n;
                if (dc * dc + dr * dr < 0.01 * n * n)
                    v = -32768.0;
                data[static_cast<size_t>(r) * n + c] = static_cast<float>(v);
            }
        RasterOptions opt;
        opt.nodata = -32768.0;
        WriteRaster(overlay, n, n, gt, &proj, data, opt);
    }
}

struct Grid {
    const char* name;
    double lon0, lat0, step;
    int n;
    bool fine;  // small enough to sample at zoom 15
};
inline std::vector<Grid> Grids() {
    return {{"site_fine", kSiteLon - 0.002, kSiteLat - 0.002, 1e-4, 41, true},
            {"overlay_edge", kSiteLon + 0.29, kSiteLat - 0.04, 2e-3, 41, true},
            {"overlay_hole", kSiteLon + 0.08, kSiteLat + 0.02, 5e-4, 33, true},
            {"regional", kSiteLon - 2.0, kSiteLat - 2.0, 0.1, 41, false}};
}

struct Point {
    const char* name;
    double lon, lat;
};
inline std::vector<Point> Points() {
    return {{"site", kSiteLon, kSiteLat},
            {"east", kSiteLon + 0.1, kSiteLat + 0.05},
            {"hole", kSiteLon + 0.098, kSiteLat + 0.037},
            {"far", -120.25, -35.5},
            {"west_wrap", -179.9, 10.0}};
}

struct Tile {
    const char* name;
    double lon, lat, half;
    int level;
};
inline Tile TileAt(const char* name, int level, double lon, double lat) {
    const double size = 16.0 / (1 << level);
    return {name, (std::floor(lon / size) + 0.5) * size, (std::floor((lat + 90) / size) + 0.5) * size - 90, 0.5 * size,
            level};
}
inline std::vector<Tile> Tiles() {
    return {TileAt("root", 0, kSiteLon, kSiteLat), TileAt("l3", 3, kSiteLon, kSiteLat),
            TileAt("l6", 6, kSiteLon, kSiteLat), TileAt("l9", 9, kSiteLon + 0.1, kSiteLat + 0.03)};
}

/// Sum, index-weighted sum, first and last value of a sequence.
template <class Seq>
std::vector<double> Checksum(const Seq& v) {
    double sum = 0, weighted = 0;
    size_t i = 0;
    for (auto x : v) {
        sum += static_cast<double>(x);
        weighted += static_cast<double>(x) * (1.0 + static_cast<double>(i % 7));
        ++i;
    }
    if (i == 0)
        return {0, 0, 0, 0};
    return {sum, weighted, static_cast<double>(*v.begin()), static_cast<double>(*(v.end() - 1))};
}

}  // namespace golden
}  // namespace planet_test

#endif
