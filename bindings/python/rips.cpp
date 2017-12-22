#include <cmath>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#include <dionysus/rips.h>

#include "simplex.h"
#include "filtration.h"

template<class T>
struct ExplicitDistances
{
    using IndexType     = int;
    using DistanceType  = float;

           ExplicitDistances(const py::array& a_):
               a(a_), n(static_cast<size_t>(1 + std::sqrt(1 + 8*a.shape(0))/2))        {}

    DistanceType operator()(int u, int v) const
    {
        if (u == v)
            return 0;

        if (u < v)              // below u must be larger than v
            std::swap(u,v);

        size_t idx = n*v - v*(v+1)/2 + u - 1 - v;
        const void* xptr = a.data(idx);
        T x = *static_cast<const T*>(xptr);

        return x;
    }

    IndexType   begin() const       { return 0; };
    IndexType   end() const         { return n; }

    const py::array&    a;
    size_t              n;
};

template<class T>
struct PairwiseDistances
{
    using IndexType     = int;
    using DistanceType  = float;

           PairwiseDistances(const py::array& a_):
               a(a_), dim(a.shape(1))       {}

    // NB: squared distance; below we pass r*r into rips.generate(...)
    DistanceType operator()(int u, int v) const
    {
        DistanceType res = 0;
        for (size_t i = 0; i < dim; ++i)
        {
            const void* xptr = a.data(u, i);
            const void* yptr = a.data(v, i);

            T x = *static_cast<const T*>(xptr);
            T y = *static_cast<const T*>(yptr);

            res += (x - y)*(x - y);
        }

        return res;
    }

    IndexType   begin() const       { return 0; };
    IndexType   end() const         { return a.shape(0); }

    const py::array&    a;
    size_t              dim;
};

template<class Distances>
PyFiltration fill_rips_(py::array a, unsigned k, double r)
{
    PyFiltration filtration;

    using Rips = dionysus::Rips<Distances, PySimplex>;
    Distances distances(a);
    Rips      rips(distances);

    rips.generate(k, r, [&filtration](PySimplex&& s) { filtration.push_back(std::move(s)); });

    typename Rips::Evaluator eval(distances);
    for (const PySimplex& s : filtration)
    {
        // in general, this is very unsafe, but we are only modifying simplex
        // data, which is not used in simplex hash, so it should be Ok overall
        PySimplex& s_ = const_cast<PySimplex&>(s);
        s_.data() = eval(s);
    }
    filtration.sort(DataDimCmp());

    return filtration;
}

PyFiltration fill_rips(py::array a, unsigned k, double r)
{
    if (a.ndim() == 2)
    {
        // PairwiseDistances returns squared distances, so we use r*r
        PyFiltration f;
        if (a.dtype().is(py::dtype::of<float>()))
            f = fill_rips_<PairwiseDistances<float>>(a,k,r*r);
        else if (a.dtype().is(py::dtype::of<double>()))
            f = fill_rips_<PairwiseDistances<double>>(a,k,r*r);
        else
            throw std::runtime_error("Unknown array dtype");

        // take square roots from simplex data
        // in general, this is very unsafe, but we are only modifying simplex
        // data, which is not used in simplex hash, so it should be Ok overall
        for (const PySimplex& s : f)
        {
            PySimplex& s_ = const_cast<PySimplex&>(s);
            s_.data() = std::sqrt(s_.data());
        }

        return f;
    } else if (a.ndim() == 1)
    {
        if (a.dtype().is(py::dtype::of<float>()))
            return fill_rips_<ExplicitDistances<float>>(a,k,r);
        else if (a.dtype().is(py::dtype::of<double>()))
            return fill_rips_<ExplicitDistances<double>>(a,k,r);
        else
            throw std::runtime_error("Unknown array dtype");
    } else
        throw std::runtime_error("Unknown input dimension: can only process 1D and 2D arrays");
}

void init_rips(py::module& m)
{
    using namespace pybind11::literals;
    m.def("fill_rips",  &fill_rips,
          "data"_a, "k"_a, "r"_a,
          "returns (sorted) filtration filled with the k-skeleton of the clique complex built on the points at distance at most r from each other");
}

