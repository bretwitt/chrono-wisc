#include "chrono_planet/dem/GeoTIFFLoader.h"

#include <cstring>
#include <iostream>
#include <memory>
#include <stdexcept>

#include <gdal_priv.h>

namespace {

OGRSpatialReference datasetSRS(GDALDataset* ds) {
    OGRSpatialReference srs;
    if (const char* wkt = ds->GetProjectionRef(); wkt && *wkt) {
        if (srs.SetFromUserInput(wkt) == OGRERR_NONE && srs.Validate() == OGRERR_NONE) {
            return srs;
        }
        std::cerr << "[GeoTIFFLoader] WKT present but could not be parsed; falling back to GetSpatialRef().\n";
    }
#if GDAL_VERSION_MAJOR >= 3
    if (const OGRSpatialReference* live = ds->GetSpatialRef()) {
        std::unique_ptr<OGRSpatialReference> copy(live->Clone());
        if (copy && copy->Validate() == OGRERR_NONE) {
            return *copy;
        }
    }
#endif
    std::cerr << "[GeoTIFFLoader] Dataset has no valid CRS, assuming EPSG:4326\n";
    srs.SetWellKnownGeogCS("WGS84");
    return srs;
}

}   // namespace

void GeoTIFFLoader::load(const std::string& filename) {
    GDALAllRegister();
    std::unique_ptr<GDALDataset, decltype(&GDALClose)> ds(
        static_cast<GDALDataset*>(GDALOpen(filename.c_str(), GA_ReadOnly)), GDALClose);
    if (!ds) {
        throw std::runtime_error("Failed to open GeoTIFF file: " + filename);
    }

    GeoTIFFLoader loaded;
    loaded.width_ = ds->GetRasterXSize();
    loaded.height_ = ds->GetRasterYSize();
    if (ds->GetGeoTransform(loaded.geoTransform_.data()) != CE_None) {
        throw std::runtime_error("Failed to get geotransform from GeoTIFF file: " + filename);
    }

    GDALRasterBand* band = ds->GetRasterBand(1);
    if (!band) {
        throw std::runtime_error("Failed to get raster band from GeoTIFF file: " + filename);
    }

    loaded.elevation_.resize(static_cast<std::size_t>(loaded.width_) * loaded.height_);
    if (band->RasterIO(GF_Read, 0, 0, loaded.width_, loaded.height_, loaded.elevation_.data(), loaded.width_, loaded.height_,
                       GDT_Float64, 0, 0) != CE_None) {
        throw std::runtime_error("RasterIO failed for GeoTIFF file: " + filename);
    }

    int hasNodata = 0;
    loaded.nodata_ = band->GetNoDataValue(&hasNodata);
    loaded.hasNodata_ = hasNodata != 0;

    loaded.srs_ = datasetSRS(ds.get());
    loaded.srs_.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    *this = std::move(loaded);
}
