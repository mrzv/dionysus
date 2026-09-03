#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "diagram.h"

#include <hera/wasserstein.h>

class WassersteinDiagramView
{
    public:
        using const_iterator = PyDiagram::const_iterator;

        explicit WassersteinDiagramView(const PyDiagram& diagram):
            diagram_(diagram)
        {}

        const_iterator begin() const      { return diagram_.begin(); }
        const_iterator end() const        { return diagram_.end(); }

    private:
        const PyDiagram& diagram_;
};

namespace hera
{
template<>
struct DiagramTraits<WassersteinDiagramView>
{
    using Container = WassersteinDiagramView;
    using PointType = PyDiagram::Point;
    using RealType  = double;

    static RealType get_x(const PointType& p)       { return p[0]; }
    static RealType get_y(const PointType& p)       { return p[1]; }
    static int      get_id(const PointType& p)      { return 0; }
};
}

using WassersteinReal = hera::DiagramTraits<WassersteinDiagramView>::RealType;

double wasserstein_distance(const PyDiagram& dgm1, const PyDiagram& dgm2, int q, double delta, double internal_p, double initial_eps, double eps_factor)
{
    hera::AuctionParams<WassersteinReal> params;
    params.wasserstein_power = q;
    params.delta = delta;
    params.internal_p = internal_p;

    if (initial_eps != 0)
        params.initial_epsilon = initial_eps;

    if (eps_factor != 0.)
        params.epsilon_common_ratio = eps_factor;

    const WassersteinDiagramView view1(dgm1);
    const WassersteinDiagramView view2(dgm2);
    return hera::wasserstein_dist(view1, view2, params);
}

void init_wasserstein_distance(py::module& m)
{
    using namespace pybind11::literals;
    m.def("wasserstein_distance",   &wasserstein_distance, "dgm1"_a, "dgm2"_a, py::arg("q") = 2,
                                                            py::arg("delta") = .01,
                                                            py::arg("internal_p") = hera::get_infinity<WassersteinReal>(),
                                                            py::arg("initial_eps") = 0.,
                                                            py::arg("eps_factor") = 0.,
          "compute Wasserstein distance between two persistence diagrams");
}
