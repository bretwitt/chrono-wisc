#ifndef QTPLANET_GEOTIFFLOADER_H
#define QTPLANET_GEOTIFFLOADER_H

#include "chrono_planet/ChApiPlanet.h"

#include <array>
#include <string>
#include <vector>

#include <ogr_spatialref.h>

// One raster band of a GeoTIFF, read whole into memory as doubles.
class CH_PLANET_API GeoTIFFLoader {
public:
    // Loads the first band. Throws std::runtime_error on failure, leaving the previous raster intact.
    void load(const std::string& filename);

    int width() const { return width_; }
    int height() const { return height_; }
    // Row-major samples, height() rows of width().
    const std::vector<double>& elevation() const { return elevation_; }
    // GDAL affine transform from (pixel, line) to dataset coordinates.
    const std::array<double, 6>& geoTransform() const { return geoTransform_; }
    const OGRSpatialReference& srs() const { return srs_; }
    // The band's nodata sentinel, if it declares one.
    bool hasNodata() const { return hasNodata_; }
    double nodata() const { return nodata_; }

private:
    int width_ = 0, height_ = 0;
    std::vector<double> elevation_;
    std::array<double, 6> geoTransform_{};
    OGRSpatialReference srs_;
    bool hasNodata_ = false;
    double nodata_ = 0.0;
};

#endif   // QTPLANET_GEOTIFFLOADER_H
