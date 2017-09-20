#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "diagram.h"

#include <../ext/hera/bottleneck/include/bottleneck.h>

double bottleneck_distance(const PyDiagram& dgm1, const PyDiagram& dgm2, double delta)
{
    if (delta == 0.0)
        return hera::bottleneckDistExact(dgm1, dgm2);
    else
        return hera::bottleneckDistApprox(dgm1, dgm2, delta);
}

void init_bottleneck_distance(py::module& m)
{
    using namespace pybind11::literals;
    m.def("bottleneck_distance",   &bottleneck_distance, "dgm1"_a, "dgm2"_a, py::arg("delta") = 0.01,
          "compute bottleneck distance between two persistence diagrams");
}
