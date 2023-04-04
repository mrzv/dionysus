#include <boost/range/adaptors.hpp>
namespace ba = boost::adaptors;

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "filtration.h"
#include "persistence.h"

PyReducedMatrix boundary(const PyFiltration& f)
{
    short prime = 3;
    PyReducedMatrix m(prime);
    m.resize(f.size());

    using Cell = PyFiltration::Cell;
    using CellChainEntry = dionysus::ChainEntry<PyReducedMatrix::Field, Cell>;
    using Entry = PyReducedMatrix::Entry;

    PyReducedMatrix::Index i = 0;
    for(auto& c : f)
    {
        m.set(i++, c.boundary(m.field()) | ba::transformed([&f,prime](const CellChainEntry& e)
                                           {
                                             short ee = e.element();
                                             if (ee > prime / 2)
                                                ee -= prime;
                                             return Entry(ee, f.index(e.index()));
                                           }));
    }

    return m;
}

PyReducedMatrix coboundary(const PyFiltration& f)
{
    short prime = 3;
    PyReducedMatrix m(prime);
    size_t n = f.size();
    m.resize(n);

    using Cell = PyFiltration::Cell;
    using CellChainEntry = dionysus::ChainEntry<PyReducedMatrix::Field, Cell>;
    using Entry = PyReducedMatrix::Entry;

    PyReducedMatrix::Index i = 0;
    for(auto& c : f)
    {
        for (auto x : c.boundary(m.field()) |
                        ba::transformed([&f,prime](const CellChainEntry& e)
                        {
                          short ee = e.element();
                          if (ee > prime / 2)
                             ee -= prime;
                          return Entry(ee, f.index(e.index()));
                        }))
        {
            m.column(n - 1 - x.index()).emplace_back(Entry { x.element(), n - 1 - i });
        }
        ++i;
    }

    for (PyReducedMatrix::Index i = 0; i < m.size(); ++i)
        m.sort(m.column(i));

    return m;
}


void init_boundary(py::module& m)
{
    m.def("boundary", &boundary, "compute boundary matrix of the filtration");
    m.def("coboundary", &coboundary, "compute coboundary matrix of the filtration");
}
