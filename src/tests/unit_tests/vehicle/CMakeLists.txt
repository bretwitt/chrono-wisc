set(TESTS
    utest_VEH_destructors
    utest_VEH_tire_contact
)

if(CH_ENABLE_MODULE_PLANET)
  set(TESTS ${TESTS} utest_VEH_PlanetTerrain)
endif()

set(UTILS
   ../ut_utils.h)

list(APPEND LIBS Chrono_core)
list(APPEND LIBS Chrono_vehicle)

message(STATUS "Add unit test programs for VEHICLE module")
build_utests(YES "${TESTS}" "vehicle" "${LIBS}" "${UTILS}")
