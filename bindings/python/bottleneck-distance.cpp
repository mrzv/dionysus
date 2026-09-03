#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "diagram.h"

#include <hera/bottleneck.h>
#include <cmath>

py::object bottleneck_distance(const PyDiagram& dgm1,
                               const PyDiagram& dgm2,
                               double delta,
                               bool compute_longest_edge)
{
    using Real = PyDiagram::Value;
    double distance;

    if (delta == 0.0)
        distance = hera::bottleneckDistExact(dgm1, dgm2, 14);
    else
        distance = hera::bottleneckDistApprox(dgm1, dgm2, delta);

    if (!compute_longest_edge)
        return py::float_(distance);

    if (distance == 0.0 || std::isinf(distance))
        return py::make_tuple(distance, py::make_tuple(-1, -1));

    hera::bt::DiagramPointSet<Real> all_a(dgm1);
    hera::bt::DiagramPointSet<Real> all_b(dgm2);
    const auto infinity_cost_edge = hera::bt::getInfinityCost(all_a, all_b, true);

    hera::bt::DiagramPointSet<Real> finite_a;
    hera::bt::DiagramPointSet<Real> finite_b;
    for (const auto& point : all_a)
        if (!point.is_infinity())
            finite_a.insert(point);
    for (const auto& point : all_b)
        if (!point.is_infinity())
            finite_b.insert(point);

    hera::bt::MatchingEdge<Real> finite_edge;
    if (delta == 0.0)
        hera::bt::bottleneckDistExact(finite_a, finite_b, 14, finite_edge, true);
    else
        hera::bt::bottleneckDistApprox(
            finite_a, finite_b, static_cast<Real>(delta), finite_edge, true);

    const auto edge_cost = [](const auto& edge)
    {
        if (edge.first.is_diagonal() && edge.second.is_normal())
            return edge.second.persistence_lp(hera::get_infinity());
        if (edge.first.is_normal() && edge.second.is_diagonal())
            return edge.first.persistence_lp(hera::get_infinity());
        return hera::bt::dist_l_inf_slow(edge.first, edge.second);
    };
    const auto finite_edge_cost = edge_cost(finite_edge);
    const auto& longest_edge = finite_edge_cost > infinity_cost_edge.cost
                                   ? finite_edge
                                   : infinity_cost_edge.edge;
    const auto index = [](const auto& point)
    {
        return point.is_diagonal() ? -1 : point.user_tag;
    };

    return py::make_tuple(distance,
                          py::make_tuple(index(longest_edge.first), index(longest_edge.second)));
}

void init_bottleneck_distance(py::module& m)
{
    using namespace pybind11::literals;
    m.def("bottleneck_distance",
          &bottleneck_distance,
          "dgm1"_a,
          "dgm2"_a,
          py::arg("delta") = 0.01,
          py::arg("compute_longest_edge") = false,
          "Compute bottleneck distance between two persistence diagrams.\n\n"
          "When compute_longest_edge is true, return "
          "(distance, (index_in_dgm1, index_in_dgm2)). A -1 index denotes "
          "the diagonal; zero and infinite distances return (-1, -1).");
}
