project                     (Dionysus)
cmake_minimum_required      (VERSION 2.8)

# Default to Release
if                          (NOT CMAKE_BUILD_TYPE)
    set                     (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    set_property            (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif                       (NOT CMAKE_BUILD_TYPE)

option                      (trace                  "Build Dionysus with trace logging"             OFF)
option                      (counters               "Build Dionysus with counters"                  OFF)
option                      (debug_zigzag           "Turn on debug routines for zigzags"            OFF)
option                      (build_examples         "Build examples"                                ON)
option                      (build_python_bindings  "Build Python bindings"                         ON)
mark_as_advanced            (debug_zigzag)

find_package                (Boost)

# Debugging
if                          (${CMAKE_BUILD_TYPE} STREQUAL "Debug" OR
                             ${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
    add_definitions         (-DDEBUG)
endif                       ()
set                         (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if                          (counters)
    add_definitions         (-DCOUNTERS)
endif                       (counters)

# Logging
if                          (trace)
    add_definitions         (-DTRACE)
endif                       (trace)

if                          (debug_zigzag)
    add_definitions         (-DDIONYSUS_ZIGZAG_DEBUG)
endif                       (debug_zigzag)

# Set includes
set                         (CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")
include_directories         (${CMAKE_CURRENT_BINARY_DIR}
                             ${CMAKE_CURRENT_SOURCE_DIR}
                             ${CMAKE_CURRENT_SOURCE_DIR}/include
                             SYSTEM ${Boost_INCLUDE_DIR})

# backward.cpp
if                          (debug)
    find_library            (LIBDW_LIBRARY NAMES dw)
    if                      (LIBDW_LIBRARY)
        set                 (DEBUG_SOURCES ${CMAKE_SOURCE_DIR}/src/backward.cpp)
        add_definitions     (-DBACKWARD_HAS_DW=1)
        set                 (libraries ${libraries} ${LIBDW_LIBRARY})
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

