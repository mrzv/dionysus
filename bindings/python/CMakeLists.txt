add_subdirectory            (pybind11)


if                          (CMAKE_LIBRARY_OUTPUT_DIRECTORY)        # to work with pip
    set                     (MODULE_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
else                        ()
    set                     (MODULE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif                       ()

file                        (MAKE_DIRECTORY "${MODULE_OUTPUT_DIRECTORY}/dionysus")
file                        (GLOB DIONYSUS_PYTHON "${CMAKE_CURRENT_SOURCE_DIR}/dionysus/*.py")
add_custom_target           (dionysus ALL ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/dionysus ${MODULE_OUTPUT_DIRECTORY}/dionysus DEPENDS ${DIONYSUS_PYTHON})

# needed to make ext/hera/wasserstein compile
include_directories         (${PROJECT_SOURCE_DIR}/ext/hera/wasserstein/include)

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
                                       bottleneck-distance.cpp
                                       wasserstein-distance.cpp)
set_target_properties       (_dionysus PROPERTIES OUTPUT_NAME dionysus/_dionysus)
