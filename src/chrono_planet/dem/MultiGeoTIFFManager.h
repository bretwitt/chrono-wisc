#ifndef QTPLANET_MULTIGEOTIFFMANAGER_H
#define QTPLANET_MULTIGEOTIFFMANAGER_H

#include "chrono_planet/ChApiPlanet.h"

#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

#include <ogr_spatialref.h>

#include "chrono_planet/dem/GeoTIFFLoader.h"

// A stack of DEMs, each valid over a range of quadtree zooms, blended so height is a pure function of
// (lon, lat, zoom) with no holes.
class CH_PLANET_API MultiGeoTIFFManager {
public:
    // Sets up the Moon 2000 lon/lat SRS projected sources are transformed to.
    MultiGeoTIFFManager();

    // Loads a GeoTIFF for zooms [min, max]. False if the file fails to load.
    [[nodiscard]] bool addSource(const std::string& filename, int minZoomInclusive, int maxZoomInclusive);

    // Height at one point, agreeing with any sampleGrid() that contains it. nullopt where no source covers it.
    [[nodiscard]] std::optional<double> sample(double lonDeg, double latDeg, int zoomLevel) const;

    struct GridSamples {
        std::vector<double> elevations;   // row-major, zero where missing
        std::vector<size_t> missing;      // indices with no valid DEM sample
    };

    // Heights over a regular lon/lat grid, row-major in j. nullopt when no source covers the grid.
    // Throws std::invalid_argument for nonpositive dimensions.
    [[nodiscard]] std::optional<GridSamples> sampleGrid(double lon0, double lat0, double stepLon, double stepLat, int nLon, int nLat,
                                                        int zoomLevel) const;

private:
    // OGR transforms must be released through GDAL, not delete.
    struct TransformDeleter {
        void operator()(OGRCoordinateTransformation* transform) const { OGRCoordinateTransformation::DestroyCT(transform); }
    };
    using TransformPtr = std::unique_ptr<OGRCoordinateTransformation, TransformDeleter>;

    // Where an overlay's data can be trusted, as a per-pixel distance to nodata or the raster edge.
    struct ValidityMask {
        std::vector<unsigned char> distancePx;   // whole pixels, saturating at 255, empty for a complete source
        float featherPx = 0.f;                   // blend ramp width in pixels
        bool complete() const { return distancePx.empty(); }
    };

    // One registered raster with its footprint, zoom range and validity mask.
    struct Source {
        std::shared_ptr<GeoTIFFLoader> loader;
        int minZoom, maxZoom;
        double minLon, maxLon, minLat, maxLat;   // degrees
        bool isGeographic;
        TransformPtr toDataset;   // lon/lat -> ds, projected sources only
        ValidityMask mask;        // complete for the base source
    };

    // Finest source valid at zoom over a rect and the coarser fallback under it. Either may be null.
    struct SelectedSources {
        const Source* fine = nullptr;
        const Source* base = nullptr;
    };
    SelectedSources selectSources(double minLon, double maxLon, double minLat, double maxLat, int zoom) const;
    // Validity mask of an overlay raster, feathered over a fixed ground distance.
    [[nodiscard]] static ValidityMask validityMaskOf(const GeoTIFFLoader& raster, double pixelMetres);
    // Blend weight in [0, 1] at a fractional pixel, from the validity mask.
    [[nodiscard]] static double weightAt(const Source& s, double colF, double rowF);
    // One source sampled onto a grid with its blend weights.
    struct WeightedGrid {
        std::vector<double> elevations;
        std::vector<double> weights;
    };
    std::optional<WeightedGrid> gridReconstruct(const Source& s, double lon0, double lat0, double stepLon, double stepLat,
                                                int nLon, int nLat) const;

    std::vector<Source> sources_;
    OGRSpatialReference moonLonLat_;
    mutable std::mutex transformMutex_;   // OGR transforms are not thread-safe
};

#endif   // QTPLANET_MULTIGEOTIFFMANAGER_H
