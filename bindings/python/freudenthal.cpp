#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#include <dionysus/rips.h>

#include "filtration.h"

struct DummyDistances
{
    using IndexType     = int;
    using DistanceType  = float;

    DistanceType operator()(int u, int v) const { return 0; }

    IndexType   begin() const       { return 0; };
    IndexType   end() const         { return 0; }
};

template<class T>
PyFiltration fill_freudenthal_(py::array a, bool reverse)
{
    PyFiltration filtration;

    using Delta        = std::vector<int>;
    using DeltaSimplex = dionysus::Simplex<Delta>;
    using Rips         = dionysus::Rips<DummyDistances, DeltaSimplex>;

    Rips::VertexContainer vertices;     // set to { 0, 0, ... }
    vertices.emplace_back();
    for (size_t i = 0; i < a.ndim(); ++i)
        vertices.back().push_back(0);

    Rips::VertexContainer candidates;
    for (size_t i = 1; i < (1 << a.ndim()); ++i)
    {
        candidates.emplace_back();
        for (size_t j = 0; j < a.ndim(); ++j)
        {
            if ((i >> j) & 1)
                candidates.back().emplace_back(1);
            else
                candidates.back().emplace_back(0);
        }
    }

    std::vector<DeltaSimplex> delta_simplices;
    Rips::bron_kerbosch(vertices, candidates, std::prev(candidates.begin()), a.ndim(),
                        [](const Delta& v1, const Delta& v2)        // neighbor test
                        {
                            bool pos = false, neg = false;
                            for (size_t i = 0; i < v1.size(); ++i)
                            {
                                auto diff = v1[i] - v2[i];
                                if (diff < 0)
                                    neg = true;
                                else if (diff > 0)
                                    pos = true;
                            }

                            return (pos ^ neg);
                        },
                        [&delta_simplices](const DeltaSimplex& s)
                        {
                            delta_simplices.push_back(s);
                        });

    // iterate over all array indices
    Delta shape(a.shape(), a.shape() + a.ndim());
    Delta v(a.ndim(), 0);
    while (v.back() <= shape.back())
    {

        for (auto& s_ : delta_simplices)
        {
            const DeltaSimplex& s = s_;
            std::vector<size_t> vertices;
            for (auto& u : s)
            {
                size_t uidx = 0;
                size_t i = 0;
                for (; i < v.size(); ++i)
                {
                    auto nbr = v[i] + u[i];
                    if (nbr >= shape[i])
                        break;
                    uidx += nbr * (a.strides(i) / sizeof(T));
                }
                if (i == v.size())
                    vertices.push_back(uidx);
            }
            if (vertices.size() == s.dimension() + 1)       // everything in bounds
                filtration.emplace_back(vertices);
        }

        size_t i = 0;
        while(i < v.size())
        {
            if (++v[i] < shape[i] || i == v.size() - 1)
                break;

            v[i] = 0;
            ++i;
        }
    }

    // fill simplex data with max/min depending on the direction
    const T* a_data = static_cast<const T*>(a.data(0));
    for (auto& s : filtration)
    {
        PySimplex& s_ = const_cast<PySimplex&>(s);
        size_t v;
        auto cmp = [a_data](size_t u, size_t v) { return *(a_data + u) < *(a_data + v); };
        if (!reverse)
            v = *std::max_element(s.begin(), s.end(), cmp);
        else
            v = *std::min_element(s.begin(), s.end(), cmp);
        s_.data() = *(a_data + v);
    }

    filtration.sort(DataDimCmp(reverse));

    return filtration;
}

PyFiltration fill_freudenthal(py::array a, bool reverse)
{
    if (a.dtype().is(py::dtype::of<float>()))
        return fill_freudenthal_<float>(a, reverse);
    else if (a.dtype().is(py::dtype::of<double>()))
        return fill_freudenthal_<double>(a, reverse);
    else
        throw std::runtime_error("Unknown array dtype");
}

void init_freudenthal(py::module& m)
{
    using namespace pybind11::literals;
    m.def("fill_freudenthal",  &fill_freudenthal,
          "data"_a, "reverse"_a = false,
          "returns (sorted) lower-star (or upper-star if ``reverse = True``) filtration filled with the Freudenthal triangulation of the grid in the array `data`");
}

