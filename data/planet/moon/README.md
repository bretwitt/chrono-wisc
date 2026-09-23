# Moon elevation data for Chrono::Planet

Elevation models used by the Moon preset (`chrono_planet/planets/moon/ChMoon.h`),
selected in code with the `moon::Dem` enum. All are LRO Lunar Orbiter Laser
Altimeter (LOLA) products (NASA, public domain), as float32 heights in meters
above the 1737.4 km reference sphere, with nodata -32768.

| `moon::Dem` | File | Resolution | Zooms | Shipped |
|---|---|---|---|---|
| `GLOBAL_LOW_RES` | `ldem_4_global.tif` | 4 px/deg (~7.6 km), global | 0-2 | yes |
| `APOLLO17_LANDING_SITE` | `ldem_1024_apollo_region.tif` | 1024 px/deg (~30 m), Apollo 17 region, equirectangular | 5-30 | yes |
| `GLOBAL` | `ldem_64_fixed.tif` | 64 px/deg (~474 m), global | 0-4 | no (480 MB) |

`ldem_4_global.tif` is `ldem_64_fixed.tif` averaged over 16 x 16 pixel blocks,
with its geotransform origin moved to the center of the first block so each
value sits where the module places it.

`GLOBAL` is too large for the repository. To use it, place a 64 px/deg LOLA
global GeoTIFF named `ldem_64_fixed.tif` in this directory (or the build tree's
copy of it). Any other DEM can be used directly through `ChGeoTiffSource`.
