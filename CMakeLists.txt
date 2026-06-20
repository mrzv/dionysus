cmake_minimum_required      (VERSION 3.20)
project                     (Dionysus VERSION 2.2.3)

include                     (CMakePackageConfigHelpers)
include                     (CTest)
include                     (GNUInstallDirs)

# Default to Release
if                          (NOT CMAKE_BUILD_TYPE)
    set                     (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    set_property            (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif                       (NOT CMAKE_BUILD_TYPE)

option                      (trace                  "Build Dionysus with trace logging"             OFF)
option                      (counters               "Build Dionysus with counters"                  OFF)
option                      (debug                  "Build Dionysus with backward.cpp debugging"    OFF)
option                      (debug_zigzag           "Turn on debug routines for zigzags"            OFF)
option                      (build_examples         "Build examples"                                ON)
option                      (build_python_bindings  "Build Python bindings"                         ON)
mark_as_advanced            (debug debug_zigzag)

find_package                (Boost CONFIG)

add_library                 (dionysus INTERFACE)
target_include_directories  (dionysus INTERFACE
                             $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                             $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_features     (dionysus INTERFACE cxx_std_14)

# Debugging
if                          (${CMAKE_BUILD_TYPE} STREQUAL "Debug" OR
                             ${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
    target_compile_definitions(dionysus INTERFACE DEBUG)
endif                       ()

if                          (counters)
    target_compile_definitions(dionysus INTERFACE COUNTERS)
endif                       (counters)

# Logging
if                          (trace)
    target_compile_definitions(dionysus INTERFACE TRACE)
endif                       (trace)

if                          (debug_zigzag)
    target_compile_definitions(dionysus INTERFACE DIONYSUS_ZIGZAG_DEBUG)
endif                       (debug_zigzag)

if                          (TARGET Boost::headers)
    target_link_libraries   (dionysus INTERFACE Boost::headers)
elseif                      (TARGET Boost::boost)
    target_link_libraries   (dionysus INTERFACE Boost::boost)
elseif                      (Boost_INCLUDE_DIRS)
    target_include_directories(dionysus SYSTEM INTERFACE ${Boost_INCLUDE_DIRS})
endif                       ()

install                     (TARGETS dionysus EXPORT DionysusTargets)
install                     (DIRECTORY include/dionysus DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install                     (EXPORT DionysusTargets
                             NAMESPACE Dionysus::
                             DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Dionysus)

configure_package_config_file(
                             cmake/DionysusConfig.cmake.in
                             ${CMAKE_CURRENT_BINARY_DIR}/DionysusConfig.cmake
                             INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Dionysus)
write_basic_package_version_file(
                             ${CMAKE_CURRENT_BINARY_DIR}/DionysusConfigVersion.cmake
                             VERSION ${PROJECT_VERSION}
                             COMPATIBILITY SameMajorVersion)
install                     (FILES
                             ${CMAKE_CURRENT_BINARY_DIR}/DionysusConfig.cmake
                             ${CMAKE_CURRENT_BINARY_DIR}/DionysusConfigVersion.cmake
                             DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Dionysus)

# backward.cpp
if                          (debug)
    find_library            (LIBDW_LIBRARY NAMES dw)
    if                      (LIBDW_LIBRARY)
        set                 (DEBUG_SOURCES ${CMAKE_SOURCE_DIR}/src/backward.cpp)
        target_compile_definitions(dionysus INTERFACE BACKWARD_HAS_DW=1)
        target_link_libraries(dionysus INTERFACE ${LIBDW_LIBRARY})
    else                    (LIBDW_LIBRARY)
        message             (STATUS "LibDW not found; backward.cpp won't be used")
    endif                   (LIBDW_LIBRARY)
endif                       (debug)

if                          (build_examples)
    add_subdirectory        (examples)
endif                       (build_examples)

if                          (build_python_bindings)
    add_subdirectory        (bindings/python)
endif                       (build_python_bindings)

if                          (BUILD_TESTING)
    add_subdirectory        (tests/cpp)
endif                       (BUILD_TESTING)
