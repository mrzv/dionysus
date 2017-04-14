#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "diagram.h"

#include <../ext/hera/wasserstein/include/wasserstein.h>

double wasserstein_distance(const PyDiagram& dgm1, const PyDiagram& dgm2, int q, double delta, double internal_p, double initial_eps, double eps_factor)
{
    return geom_ws::wassersteinDist(dgm1, dgm2, q, delta, internal_p, initial_eps, eps_factor);
}

void init_distances(py::module& m)
{
    using namespace pybind11::literals;
    m.def("wasserstein_distance",   &wasserstein_distance, "dgm1"_a, "dgm2"_a, py::arg("q") = 2,
                                                            py::arg("delta") = .01,
                                                            py::arg("internal_p") = std::numeric_limits<double>::infinity(),
                                                            py::arg("initial_eps") = 0.,
                                                            py::arg("eps_factor") = 0.,
          "compute Wasserstein distance between two persistence diagrams");
}
