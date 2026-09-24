#=============================================================================
# CMake configuration file for demos of the Chrono Planet module.
# Requires the Vehicle and VSG modules.
#=============================================================================

set(DEMOS
    demo_PLANET_Surface
    demo_PLANET_TerrainVSG
    demo_PLANET_Flyover
    demo_PLANET_Viper_Rigid
    demo_PLANET_Viper_SCM
)

if(CH_ENABLE_MODULE_VEHICLE_MODELS)
    set(DEMOS ${DEMOS}
        demo_PLANET_WheeledVehicle
    )
endif()

set(UTILS "PlanetBoulderField.h" "PlanetDemoSetup.h" "PlanetHazardAvoidance.h")

#--------------------------------------------------------------

list(APPEND LIBS Chrono_core)
list(APPEND LIBS Chrono_planet)
list(APPEND LIBS Chrono_planet_vsg)
list(APPEND LIBS Chrono_vehicle)
list(APPEND LIBS Chrono_vehicle_vsg)
list(APPEND LIBS ChronoModels_robot)
list(APPEND LIBS Chrono_vsg)

if(CH_ENABLE_MODULE_VEHICLE_MODELS)
    list(APPEND LIBS ChronoModels_vehicle)
endif()

#--------------------------------------------------------------

message(STATUS "Add demo programs for PLANET module")

foreach(PROGRAM ${DEMOS})
    add_executable(${PROGRAM} ${PROGRAM}.cpp ${UTILS})
    source_group("" FILES ${PROGRAM}.cpp ${UTILS})

    if(APPLE)
      set_target_properties(${PROGRAM} PROPERTIES MACOSX_BUNDLE_INFO_PLIST ${MACOSX_BUNDLE_INFO_PLIST})
    endif()
    if(MSVC)
      set_target_properties(${PROGRAM} PROPERTIES MSVC_RUNTIME_LIBRARY ${CH_MSVC_RUNTIME_LIBRARY})
    endif()
    set_property(TARGET ${PROGRAM} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROGRAM}>")
    target_link_libraries(${PROGRAM} ${LIBS})

    install(TARGETS ${PROGRAM} CONFIGURATIONS Release DESTINATION ${CH_INSTALL_DEMO})
endforeach()

# Chrono::Sensor demos over the planet terrain, shaded with the Hapke BRDF.
if(CH_ENABLE_MODULE_SENSOR)
    foreach(PROGRAM demo_PLANET_SensorHapke demo_PLANET_Viper_SCM_Sensor demo_PLANET_MoonDescent)
        add_executable(${PROGRAM} ${PROGRAM}.cpp ${UTILS})
        source_group("" FILES ${PROGRAM}.cpp ${UTILS})
        if(MSVC)
          set_target_properties(${PROGRAM} PROPERTIES MSVC_RUNTIME_LIBRARY ${CH_MSVC_RUNTIME_LIBRARY})
        endif()
        set_property(TARGET ${PROGRAM} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROGRAM}>")
        target_link_libraries(${PROGRAM} Chrono_core Chrono_planet Chrono_vehicle ChronoModels_robot Chrono_sensor)
        install(TARGETS ${PROGRAM} CONFIGURATIONS Release DESTINATION ${CH_INSTALL_DEMO})
    endforeach()
endif()
