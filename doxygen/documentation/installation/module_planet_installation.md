Install the PLANET module {#module_planet_installation}
==========================

[TOC]

Chrono::Planet is an optional module that provides planetary terrain from
GeoTIFF digital elevation models: a stack of DEMs sampled at any longitude,
latitude and level of detail, procedural crater, rock and roughness fields
below the DEM resolution, and a quadtree of terrain tiles that streams around
a moving viewpoint. Chrono::Vehicle builds two terrain models on it:
`PlanetTerrain` (a rigid collision patch that follows the vehicle) and
`PlanetSCMTerrain` (the SCM soil model over the unbounded procedural surface).


## Features

- GeoTIFF loading through GDAL, with several DEMs blended by zoom range and
  feathered where one overlays another.
- A deterministic surface (`chrono::planet::ChPlanetSurface`) for physics,
  evaluated at a fixed zoom so wheels and rendering share one height field.
- A local east-north-up frame at a landing site (`chrono::planet::ChSiteFrame`).
- Procedural craters, boulders and micro-roughness that are pure functions of
  longitude and latitude, so every consumer agrees.
- A quadtree of tiles (`QuadtreeWorld`) with a prefetch worker and a tile mesh
  builder, ready for a rendering back end.


## Requirements

- To **build** this module you need:
    - the [GDAL](https://gdal.org) library and headers (any 3.x)
    - optionally OpenMP, which the grid loops use when Chrono is built with it

- To **run** applications based on this module you need:
    - the GDAL shared library
    - one or more GeoTIFF elevation models; with none, the surface is a
      procedural fallback plus craters, which is enough for the demos


## Building instructions

1. Install GDAL (on Ubuntu, `libgdal-dev`).

2. Repeat the instructions for the [full installation](@ref tutorial_install_chrono)

3. During CMake configuration, set `CH_ENABLE_MODULE_PLANET` to 'on'.
   CMake's own FindGDAL locates the library; set `GDAL_DIR` or `GDAL_ROOT`
   if it is not in a default location.

4. Optional settings:
    - `CH_PLANET_SIMD` selects the instruction set for the terrain grid loops
      (`avx2`, the default, `native`, or `off`). Contraction is always off so
      every build produces the same heights.
    - `CH_PLANET_ASSET_ROOT` bakes in a directory holding a `resources/` tree,
      so DEM paths given relative to it resolve regardless of the working
      directory.

5. Press 'Generate' and build as usual.

To use the vehicle terrains, also enable `CH_ENABLE_MODULE_VEHICLE`.

To draw the quadtree terrain, also enable `CH_ENABLE_MODULE_VSG`: this adds
the `Chrono_planet_vsg` library with `ChPlanetVisualizationVSG`, a plugin for
`ChVisualSystemVSG` that streams the tiles around the camera, and the
`demo_PLANET_TerrainVSG` demo.


## How to use it

- Look at the [API section](group__planet__module.html) of this module for documentation about classes and functions.
- Look at the C++ source of [demos](@ref tutorial_root) to learn how to use the functions of this module;
  `demo_VEH_PlanetTerrain` runs both terrain models headless with no data files.
