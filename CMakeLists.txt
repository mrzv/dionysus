project                     (Dionysus)
cmake_minimum_required      (VERSION 2.8)

option                      (debug              "Build Dionysus with debugging"                 OFF)
option                      (optimize           "Build Dionysus with optimization"              ON)
option                      (trace              "Build Dionysus with trace logging"             OFF)
option                      (counters           "Build Dionysus with counters"                  OFF)
option                      (debug_zigzag       "Turn on debug routines for zigzags"            OFF)
mark_as_advanced            (debug_zigzag)

# Debugging
if                          (debug)
    if                      (optimize)
            set             (cxx_flags                  ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
    else                    (optimize)
            set             (cxx_flags                  ${CMAKE_CXX_FLAGS_DEBUG})
    endif                   (optimize)
    add_definitions         (-DDEBUG)
else                        (debug)
    if                      (optimize)
            set             (cxx_flags                  ${CMAKE_CXX_FLAGS_RELEASE})
    else                    (optimize)
            set             (cxx_flags                  ${CMAKE_CXX_FLAGS})
    endif                   (optimize)
endif                       (debug)
add_definitions             (${cxx_flags} -std=c++11)

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

add_subdirectory            (examples)


