set(TESTS
    utest_PLANET_Generic
    utest_PLANET_MoonRegression
)

set(UTILS
    planet_test_utils.h)

# The tests write their own GeoTIFFs, so they use GDAL directly.
find_package(GDAL REQUIRED)

list(APPEND LIBS Chrono_core)
list(APPEND LIBS Chrono_planet)
list(APPEND LIBS GDAL::GDAL)

message(STATUS "Add unit test programs for PLANET module")
build_utests(YES "${TESTS}" "planet" "${LIBS}" "${UTILS}")
