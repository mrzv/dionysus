add_subdirectory            (pybind11)


if                          (CMAKE_LIBRARY_OUTPUT_DIRECTORY)        # to work with pip
    set                     (MODULE_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
else                        ()
    set                     (MODULE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif                       ()

file                        (MAKE_DIRECTORY "${MODULE_OUTPUT_DIRECTORY}/dionysus")
file                        (GLOB DIONYSUS_PYTHON "${CMAKE_CURRENT_SOURCE_DIR}/dionysus/*.py")
add_custom_target           (dionysus ALL ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/dionysus ${MODULE_OUTPUT_DIRECTORY}/dionysus DEPENDS ${DIONYSUS_PYTHON})

pybind11_add_module         (_dionysus dionysus.cpp
                                       filtration.cpp
                                       simplex.cpp
                                       field.cpp
                                       rips.cpp
                                       freudenthal.cpp
                                       persistence.cpp
                                       diagram.cpp
                                       omni-field-persistence.cpp
                                       cohomology-persistence.cpp
                                       zigzag-persistence.cpp
                                       distances.cpp)
target_link_libraries       (_dionysus PRIVATE wasserstein)     # CMake will find the target in hera/wasserstein;
                                                                # PRIVATE to keep consistent with target_link_libraries in pybind11_add_module
set_target_properties       (_dionysus PROPERTIES OUTPUT_NAME dionysus/_dionysus)
