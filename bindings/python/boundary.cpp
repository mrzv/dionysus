#include <boost/range/adaptors.hpp>
namespace ba = boost::adaptors;

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "filtration.h"
#include "persistence.h"

template<class PyFiltration>
PyMatrixFiltration boundary(const PyFiltration& f)
{
    short prime = 3;
    PyReducedMatrix m(prime);
    Dimensions dimensions;
    Values values;
    m.resize(f.size());
    dimensions.resize(f.size());
    values.resize(f.size());

    using Cell = typename PyFiltration::Cell;
    using CellChainEntry = dionysus::ChainEntry<PyReducedMatrix::Field, Cell>;
    using Entry = PyReducedMatrix::Entry;

    PyReducedMatrix::Index i = 0;
    for(auto& c : f)
    {
        dimensions[i] = c.dimension();
        values[i] = c.data();
        m.set(i, c.boundary(m.field()) | ba::transformed([&f,prime,i](const CellChainEntry& e)
                                         {
                                           short ee = e.element();
                                           if (ee > prime / 2)
                                              ee -= prime;
                                           return Entry(ee, f.index(e.index(),i));
                                         }));
        ++i;
    }

    return PyMatrixFiltration(std::move(m),dimensions,values);
}

template<class PyFiltration>
PyMatrixFiltration coboundary(const PyFiltration& f)
{
    short prime = 3;
    PyReducedMatrix m(prime);
    Dimensions dimensions;
    Values values;

    using Cell = typename PyFiltration::Cell;
    using CellChainEntry = dionysus::ChainEntry<PyReducedMatrix::Field, Cell>;
    using Entry = PyReducedMatrix::Entry;
    using Index = PyReducedMatrix::Index;

    Index n = static_cast<Index>(f.size());
    m.resize(n);
    dimensions.resize(n);
    values.resize(n);

    Index i = 0;
    for(auto& c : f)
    {
        dimensions[n - 1 - i] = c.dimension();
        values[n - 1 - i] = c.data();
        for (auto x : c.boundary(m.field()) |
                        ba::transformed([&f,prime,i](const CellChainEntry& e)
                        {
                          short ee = e.element();
                          if (ee > prime / 2)
                             ee -= prime;
                          return Entry(ee, f.index(e.index(),i));
                        }))
        {
            m.column(n - 1 - x.index()).emplace_back(Entry { x.element(), n - 1 - i });
        }
        ++i;
    }

    for (PyReducedMatrix::Index i = 0; i < m.size(); ++i)
        m.sort(m.column(i));

    return PyMatrixFiltration(std::move(m),dimensions,values);
}


void init_boundary(py::module& m)
{
    m.def("boundary", &boundary<PyFiltration>, "compute boundary matrix of the filtration");
    m.def("coboundary", &coboundary<PyFiltration>, "compute coboundary matrix of the filtration");
    m.def("boundary", &boundary<PyMultiFiltration>, "compute boundary matrix of the filtration");
    m.def("coboundary", &coboundary<PyMultiFiltration>, "compute coboundary matrix of the filtration");
}
