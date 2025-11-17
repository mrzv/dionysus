set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

file                        (GLOB DIONYSUS_PYTHON "${CMAKE_CURRENT_SOURCE_DIR}/dionysus/*.py")
install                     (FILES ${DIONYSUS_PYTHON} DESTINATION dionysus)

# needed to make ext/hera compile
include_directories         (${PROJECT_SOURCE_DIR}/ext/hera/include)

pybind11_add_module         (_dionysus dionysus.cpp
                                       filtration.cpp
                                       simplex.cpp
                                       field.cpp
                                       rips.cpp
                                       freudenthal.cpp
                                       persistence.cpp
                                       boundary.cpp
                                       diagram.cpp
                                       omni-field-persistence.cpp
                                       cohomology-persistence.cpp
                                       zigzag-persistence.cpp
                                       bottleneck-distance.cpp
                                       wasserstein-distance.cpp)

install                     (TARGETS _dionysus DESTINATION dionysus)
