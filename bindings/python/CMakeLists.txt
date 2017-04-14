list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
find_package(PythonLibsNew ${PYBIND11_PYTHON_VERSION} REQUIRED)

include(CheckCXXCompilerFlag)

if(NOT MSVC AND NOT PYBIND11_CPP_STANDARD)
  check_cxx_compiler_flag("-std=c++14" HAS_CPP14_FLAG)
  check_cxx_compiler_flag("-std=c++11" HAS_CPP11_FLAG)

  if (HAS_CPP14_FLAG)
    set(PYBIND11_CPP_STANDARD -std=c++14)
  elseif (HAS_CPP11_FLAG)
    set(PYBIND11_CPP_STANDARD -std=c++11)
  else()
    message(FATAL_ERROR "Unsupported compiler -- pybind11 requires C++11 support!")
  endif()

  set(PYBIND11_CPP_STANDARD ${PYBIND11_CPP_STANDARD} CACHE STRING
      "C++ standard flag, e.g. -std=c++11 or -std=c++14. Defaults to latest available.")
endif()

# Cache variables so pybind11_add_module can be used in parent projects
set(PYBIND11_INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/include" CACHE INTERNAL "")
set(PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS} CACHE INTERNAL "")
set(PYTHON_LIBRARIES ${PYTHON_LIBRARIES} CACHE INTERNAL "")
set(PYTHON_MODULE_PREFIX ${PYTHON_MODULE_PREFIX} CACHE INTERNAL "")
set(PYTHON_MODULE_EXTENSION ${PYTHON_MODULE_EXTENSION} CACHE INTERNAL "")

# Build a Python extension module:
# pybind11_add_module(<name> source1 [source2 ...])
#
function(pybind11_add_module target_name)
  add_library(${target_name} MODULE ${ARGN})
  target_include_directories(${target_name} PUBLIC ${PYBIND11_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS})

  # The prefix and extension are provided by FindPythonLibsNew.cmake
  set_target_properties(${target_name} PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}")
  set_target_properties(${target_name} PROPERTIES SUFFIX "${PYTHON_MODULE_EXTENSION}")

  if(WIN32 OR CYGWIN)
    # Link against the Python shared library on Windows
    target_link_libraries(${target_name} PRIVATE ${PYTHON_LIBRARIES})
  elseif(APPLE)
    # It's quite common to have multiple copies of the same Python version
    # installed on one's system. E.g.: one copy from the OS and another copy
    # that's statically linked into an application like Blender or Maya.
    # If we link our plugin library against the OS Python here and import it
    # into Blender or Maya later on, this will cause segfaults when multiple
    # conflicting Python instances are active at the same time (even when they
    # are of the same version).

    # Windows is not affected by this issue since it handles DLL imports
    # differently. The solution for Linux and Mac OS is simple: we just don't
    # link against the Python library. The resulting shared library will have
    # missing symbols, but that's perfectly fine -- they will be resolved at
    # import time.

    target_link_libraries(${target_name} PRIVATE "-undefined dynamic_lookup")
  endif()

  if(NOT MSVC)
    # Make sure C++11/14 are enabled
    target_compile_options(${target_name} PUBLIC ${PYBIND11_CPP_STANDARD})

    # Enable link time optimization and set the default symbol
    # visibility to hidden (very important to obtain small binaries)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)
    if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
      # Check for Link Time Optimization support (GCC/Clang)
      check_cxx_compiler_flag("-flto" HAS_LTO_FLAG)
      if(HAS_LTO_FLAG AND NOT CYGWIN)
        target_compile_options(${target_name} PRIVATE -flto)
      endif()

      # Intel equivalent to LTO is called IPO
      if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        check_cxx_compiler_flag("-ipo" HAS_IPO_FLAG)
        if(HAS_IPO_FLAG)
          target_compile_options(${target_name} PRIVATE -ipo)
        endif()
      endif()

      # Default symbol visibility
      target_compile_options(${target_name} PRIVATE "-fvisibility=hidden")

      # Strip unnecessary sections of the binary on Linux/Mac OS
      if(CMAKE_STRIP)
        if(APPLE)
          add_custom_command(TARGET ${target_name} POST_BUILD
                             COMMAND ${CMAKE_STRIP} -u -r $<TARGET_FILE:${target_name}>)
        else()
          add_custom_command(TARGET ${target_name} POST_BUILD
                             COMMAND ${CMAKE_STRIP} $<TARGET_FILE:${target_name}>)
        endif()
      endif()
    endif()
  elseif(MSVC)
    # /MP enables multithreaded builds (relevant when there are many files), /bigobj is
    # needed for bigger binding projects due to the limit to 64k addressable sections
    target_compile_options(${target_name} PRIVATE /MP /bigobj)

    # Enforce link time code generation on MSVC, except in debug mode
    target_compile_options(${target_name} PRIVATE $<$<NOT:$<CONFIG:Debug>>:/GL>)

    # Fancy generator expressions don't work with linker flags, for reasons unknown
    set_property(TARGET ${target_name} APPEND_STRING PROPERTY LINK_FLAGS_RELEASE /LTCG)
    set_property(TARGET ${target_name} APPEND_STRING PROPERTY LINK_FLAGS_MINSIZEREL /LTCG)
    set_property(TARGET ${target_name} APPEND_STRING PROPERTY LINK_FLAGS_RELWITHDEBINFO /LTCG)
  endif()
endfunction()

file                        (MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/dionysus")
file                        (GLOB DIONYSUS_PYTHON "${CMAKE_CURRENT_SOURCE_DIR}/dionysus/*.py")
add_custom_target           (dionysus ALL ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/dionysus ${CMAKE_CURRENT_BINARY_DIR}/dionysus DEPENDS ${DIONYSUS_PYTHON})

pybind11_add_module         (_dionysus dionysus.cpp
                                       filtration.cpp
                                       simplex.cpp
                                       field.cpp
                                       rips.cpp
                                       persistence.cpp
                                       diagram.cpp
                                       omni-field-persistence.cpp
                                       cohomology-persistence.cpp
                                       distances.cpp)
target_link_libraries       (_dionysus PRIVATE wasserstein)     # CMake will find the target in hera/wasserstein;
                                                                # PRIVATE to keep consistent with target_link_libraries in pybind11_add_module
set_target_properties       (_dionysus PROPERTIES OUTPUT_NAME dionysus/_dionysus)
