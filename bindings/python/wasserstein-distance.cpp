#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "diagram.h"

#include <../ext/hera/wasserstein/include/wasserstein.h>

double wasserstein_distance(const PyDiagram& dgm1, const PyDiagram& dgm2, int q, double delta, double internal_p, double initial_eps, double eps_factor)
{
    hera::AuctionParams<PyDiagram::Value> params;
    params.wasserstein_power = q;
    params.delta = delta;
    params.internal_p = internal_p;

    if (initial_eps != 0)
        params.initial_epsilon = initial_eps;

    if (eps_factor != 0.)
        params.epsilon_common_ratio = eps_factor;

    return hera::wasserstein_dist(dgm1, dgm2, params);
}

void init_wasserstein_distance(py::module& m)
{
    using namespace pybind11::literals;
    m.def("wasserstein_distance",   &wasserstein_distance, "dgm1"_a, "dgm2"_a, py::arg("q") = 2,
                                                            py::arg("delta") = .01,
                                                            py::arg("internal_p") = hera::get_infinity<PyDiagram::Value>(),
                                                            py::arg("initial_eps") = 0.,
                                                            py::arg("eps_factor") = 0.,
          "compute Wasserstein distance between two persistence diagrams");
}
