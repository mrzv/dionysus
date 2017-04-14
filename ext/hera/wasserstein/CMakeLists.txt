project                     (wasserstein)
cmake_minimum_required      (VERSION 2.8.9)
include (GenerateExportHeader)

# Default to Release

if                          (NOT CMAKE_BUILD_TYPE)
    set                     (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif                       (NOT CMAKE_BUILD_TYPE)

# Boost
find_package                (Boost)
include_directories         (${CMAKE_CURRENT_SOURCE_DIR}/include
                             SYSTEM ${Boost_INCLUDE_DIR})

if(NOT WIN32)
    add_definitions(-std=c++11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 ")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -O3 -g")
endif(NOT WIN32)

file(GLOB WS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

add_library(wasserstein ${WS_SOURCES})

if (WIN32)
    GENERATE_EXPORT_HEADER(wasserstein
               BASE_NAME wasserstein
               EXPORT_MACRO_NAME wasserstein_EXPORT
               EXPORT_FILE_NAME wasserstein_export.h
               STATIC_DEFINE wasserstein_BUILT_AS_STATIC)
endif(WIN32)

add_executable(wasserstein_dist ${CMAKE_CURRENT_SOURCE_DIR}/example/wasserstein_dist.cpp)
target_link_libraries(wasserstein_dist PUBLIC wasserstein)
