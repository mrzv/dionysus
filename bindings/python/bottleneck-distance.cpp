#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "diagram.h"

#include <../ext/hera/bottleneck/include/bottleneck.h>

py::object
bottleneck_distance(const PyDiagram& dgm1, const PyDiagram& dgm2, double delta, bool compute_longest_edge = true)
{
    hera::bt::MatchingEdge<double> longest_edge;

    // dionysus uses float,
    // safer to use double in Hera code
    // so we copy input diagrams
    std::vector<std::pair<double, double>> dgm_a, dgm_b;
    py::object result;
    dgm_a.reserve(dgm1.size());
    dgm_b.reserve(dgm2.size());

    for (size_t i = 0; i < dgm1.size(); ++i) {
        dgm_a.emplace_back(dgm1[i].birth(), dgm1[i].death());
    }
    for (size_t i = 0; i < dgm2.size(); ++i) {
        dgm_b.emplace_back(dgm2[i].birth(), dgm2[i].death());
    }

    double d;

    if (delta == 0.0) {
        int dec_precision = 14;
        d = hera::bottleneckDistExact(dgm_a, dgm_b, dec_precision, longest_edge, compute_longest_edge);
    } else {
        d = hera::bottleneckDistApprox(dgm_a, dgm_b, delta, longest_edge, compute_longest_edge);
    }
    if (compute_longest_edge) {
        result = py::cast(std::make_pair(d,
                                         std::make_pair(longest_edge.first.get_user_id(),
                                                        longest_edge.second.get_user_id())));
    } else {
        result = py::cast(d);
    }
    return result;
}


void init_bottleneck_distance(py::module& m)
{
    using namespace pybind11::literals;
    m.def("bottleneck_distance", &bottleneck_distance, "dgm1"_a, "dgm2"_a, py::arg("delta") = 0.01, py::arg("compute_longest_edge") = false,
          "compute bottleneck distance between two persistence diagrams and longest edge");
}
