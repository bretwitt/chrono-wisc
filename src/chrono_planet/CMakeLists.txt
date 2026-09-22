#=============================================================================
# CMake configuration file for the Chrono Planet module
#
# Cannot be used stand-alone (it's loaded by CMake config. file in parent dir.)
#=============================================================================

option(CH_ENABLE_MODULE_PLANET "Enable the Chrono Planet module (GeoTIFF DEM terrain with quadtree LOD)" OFF)

# Return now if this module is not enabled
if(NOT CH_ENABLE_MODULE_PLANET)
  return()
endif()

message(STATUS "\n==== Chrono Planet module ====\n")

#-----------------------------------------------------------------------------
# GDAL reads the GeoTIFFs. CMake's own FindGDAL provides the GDAL::GDAL target.
#-----------------------------------------------------------------------------

find_package(GDAL REQUIRED)
message(STATUS "GDAL version: ${GDAL_VERSION}")

#-----------------------------------------------------------------------------
# ISA for the terrain core's grid loops. Contraction stays off so every build
# sums the same bits; "native" is faster but not reproducible across machines.
#-----------------------------------------------------------------------------

set(CH_PLANET_SIMD "avx2" CACHE STRING "Planet module grid-loop ISA: avx2, native or off")
set_property(CACHE CH_PLANET_SIMD PROPERTY STRINGS avx2 native off)
if(CH_PLANET_SIMD STREQUAL "native")
  set(CH_PLANET_SIMD_FLAGS "$<$<CXX_COMPILER_ID:GNU,Clang>:-march=native;-ffp-contract=off>")
elseif(CH_PLANET_SIMD STREQUAL "avx2")
  set(CH_PLANET_SIMD_FLAGS "$<$<CXX_COMPILER_ID:GNU,Clang>:-mavx2;-ffp-contract=off>")
else()
  set(CH_PLANET_SIMD_FLAGS "$<$<CXX_COMPILER_ID:GNU,Clang>:-ffp-contract=off>")
endif()

# Optional directory holding the renderer's shaders/ and resources/, baked in so
# relative DEM paths resolve without depending on the working directory.
set(CH_PLANET_ASSET_ROOT "" CACHE PATH "Directory containing the terrain shaders/ and resources/ trees (optional)")

#-----------------------------------------------------------------------------

set(Chrono_PLANET_CORE_FILES
    core/AssetPaths.h
    core/AssetPaths.cpp
    core/BenchProfiler.h
    core/BenchProfiler.cpp
    core/FilterChain.h
    core/MathUtil.h
    core/MathUtil.cpp
    core/Parallel.h
    core/Planet.h
    core/SphereMath.h
    core/SphereMath.cpp
)
source_group("core" FILES ${Chrono_PLANET_CORE_FILES})

set(Chrono_PLANET_DEM_FILES
    dem/GeoTIFFLoader.h
    dem/GeoTIFFLoader.cpp
    dem/MultiGeoTIFFManager.h
    dem/MultiGeoTIFFManager.cpp
)
source_group("dem" FILES ${Chrono_PLANET_DEM_FILES})

set(Chrono_PLANET_LOD_FILES
    lod/CoordinateSystems.h
    lod/Quadtree.h
    lod/Quadtree.cpp
    lod/QuadtreeTile.h
    lod/QuadtreeTile.cpp
    lod/QuadtreeWorld.h
    lod/QuadtreeWorld.cpp
    lod/SphericalCoordinates.h
    lod/SphericalCoordinates.cpp
    lod/TileBuildWorker.h
    lod/TileBuildWorker.cpp
    lod/TileMeshBuilder.h
    lod/TileMeshBuilder.cpp
    lod/TileMetadata.h
    lod/TileMetadata.cpp
)
source_group("lod" FILES ${Chrono_PLANET_LOD_FILES})

set(Chrono_PLANET_PROCEDURAL_FILES
    procedural/CraterField.h
    procedural/CraterField.cpp
    procedural/FieldGrid.h
    procedural/FieldHash.h
    procedural/Perlin.h
    procedural/Perlin.cpp
    procedural/RockField.h
    procedural/RockField.cpp
    procedural/RockMeshes.h
    procedural/RockMeshes.cpp
    procedural/RockSdf.h
    procedural/RockSdf.cpp
    procedural/SurfaceRoughness.h
    procedural/SurfaceRoughness.cpp
)
source_group("procedural" FILES ${Chrono_PLANET_PROCEDURAL_FILES})

set(Chrono_PLANET_SHADING_FILES
    shading/Hapke.h
    shading/Hapke.cpp
)
source_group("shading" FILES ${Chrono_PLANET_SHADING_FILES})

set(Chrono_PLANET_FILES
    ChApiPlanet.h
    ChSiteFrame.h
    ChPlanetSurface.h
    ChPlanetSurface.cpp
)
source_group("" FILES ${Chrono_PLANET_FILES})

# The grid loops in these files only vectorise if the compiler may hoist
# floating-point work out of a conditional and skip the errno check on sqrt.
# Neither changes a result and nothing here reads errno.
set_source_files_properties(
    procedural/CraterField.cpp
    procedural/RockField.cpp
    procedural/SurfaceRoughness.cpp
    dem/MultiGeoTIFFManager.cpp
    lod/SphericalCoordinates.cpp
    lod/TileMeshBuilder.cpp
    PROPERTIES COMPILE_OPTIONS "$<$<CXX_COMPILER_ID:GNU,Clang>:-fno-math-errno;-fno-trapping-math>"
)

#-----------------------------------------------------------------------------
# In most cases, you do not need to edit the lines below.

add_library(Chrono_planet
    ${Chrono_PLANET_FILES}
    ${Chrono_PLANET_CORE_FILES}
    ${Chrono_PLANET_DEM_FILES}
    ${Chrono_PLANET_LOD_FILES}
    ${Chrono_PLANET_PROCEDURAL_FILES}
    ${Chrono_PLANET_SHADING_FILES}
)
add_library(Chrono::planet ALIAS Chrono_planet)

set_target_properties(Chrono_planet PROPERTIES DEBUG_POSTFIX ${CH_DEBUG_POSTFIX})

if(CH_WHOLE_PROG_OPT)
  set_target_properties(Chrono_planet PROPERTIES COMPILE_FLAGS "/GL")
  set_target_properties(Chrono_planet PROPERTIES LINK_FLAGS "/LTCG")
endif()

if(CH_STATIC)
  set_target_properties(Chrono_planet PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

if(MSVC)
  set_target_properties(Chrono_planet PROPERTIES MSVC_RUNTIME_LIBRARY ${CH_MSVC_RUNTIME_LIBRARY})
endif()

target_compile_definitions(Chrono_planet PRIVATE $<$<COMPILE_LANGUAGE:CXX>:CH_API_COMPILE_PLANET>)
if(CH_PLANET_ASSET_ROOT)
  target_compile_definitions(Chrono_planet PRIVATE QTPLANET_SOURCE_ROOT="${CH_PLANET_ASSET_ROOT}")
endif()

target_compile_options(Chrono_planet PRIVATE ${CH_PLANET_SIMD_FLAGS})

# GDAL is public: the DEM headers expose OGRSpatialReference members.
target_link_libraries(Chrono_planet PUBLIC Chrono_core GDAL::GDAL)

install(TARGETS Chrono_planet
        EXPORT ChronoTargets
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include/chrono_planet)

#-------------------------------------------------------------------------------
# Install files
#-------------------------------------------------------------------------------

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
        DESTINATION include/chrono_planet
        FILES_MATCHING PATTERN "*.h")

# Dependencies on other Chrono components, for the exported ChronoConfig.
set(DEPENDENCIES_PLANET "" PARENT_SCOPE)
