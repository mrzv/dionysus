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

std::pair<double, std::pair<int, int>> bottleneck_distance_with_edge(const PyDiagram& dgm1, const PyDiagram& dgm2, double delta)
{
    hera::bt::MatchingEdge<float> longest_edge;
    bool compute_longest_edge = true;

    if (delta == 0.0) {
        int dec_precision = 14;
        float d = hera::bottleneckDistExact(dgm1, dgm2, dec_precision, longest_edge, compute_longest_edge);

        return std::make_pair(d,
                              std::make_pair(longest_edge.first.get_user_id(),
                                  longest_edge.second.get_user_id()));
    } else {
        float d = hera::bottleneckDistApprox(dgm1, dgm2, delta, longest_edge, compute_longest_edge);

        return std::make_pair(d,
                              std::make_pair(longest_edge.first.get_user_id(),
                                  longest_edge.second.get_user_id()));
    }
}


void init_bottleneck_distance(py::module& m)
{
    using namespace pybind11::literals;
    m.def("bottleneck_distance",   &bottleneck_distance, "dgm1"_a, "dgm2"_a, py::arg("delta") = 0.01,
          "compute bottleneck distance between two persistence diagrams");
    m.def("bottleneck_distance_with_edge",   &bottleneck_distance_with_edge, "dgm1"_a, "dgm2"_a, py::arg("delta") = 0.01,
          "compute bottleneck distance between two persistence diagrams and the longest edge");
}
