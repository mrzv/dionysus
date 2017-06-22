#include <boost/range/adaptors.hpp>
namespace ba = boost::adaptors;

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "filtration.h"
#include "diagram.h"
#include "omni-field-persistence.h"

PyOmniFieldPersistence
omnifield_homology_persistence(const PyFiltration& filtration)
{
    PyOmniFieldPersistence persistence;
    for(auto& s : filtration)
    {
        using SimplexChainEntry = dionysus::ChainEntry<PyOmniFieldPersistence::Field, PySimplex>;
        using ChainEntry        = dionysus::ChainEntry<PyOmniFieldPersistence::Field, PyOmniFieldPersistence::Index>;
        persistence.add(s.boundary(persistence.field()) |
                                                 ba::transformed([&filtration](const SimplexChainEntry& e)
                                                 { return ChainEntry(e.element(), filtration.index(e.index())); }));
    }
    return persistence;
}

std::vector<PyDiagram>
py_init_omni_diagrams(const PyOmniFieldPersistence& persistence, const PyFiltration& f, PyOmniFieldPersistence::BaseElement p)
{
    return init_diagrams(prime_adapter(persistence, p), f,
                         [](const PySimplex& s)                         { return s.data(); },        // value
                         [](PyOmniFieldPersistence::Index i) -> PyIndex { return i; });              // data
}

PYBIND11_MAKE_OPAQUE(PyOmniFieldPersistence::ZpChain);      // persistence.cpp provides a binding for Chain, which is exactly what this is

void init_omnifield_persistence(py::module& m)
{
    using namespace pybind11::literals;
    m.def("omnifield_homology_persistence",   &omnifield_homology_persistence, "filtration"_a,
          "compute homology persistence of the filtration (pair simplices) over all fields at once");

    m.def("init_diagrams",      &py_init_omni_diagrams,  "ofp"_a, "f"_a, "p"_a,  "initialize diagrams for a specific prime from omnifield persistence and filtration");

    using Index         = PyOmniFieldPersistence::Index;
    using BaseElement   = PyOmniFieldPersistence::BaseElement;
    py::class_<PyOmniFieldPersistence>(m, "OmniFieldPersistence", "compact composition of multiple reduced matrices")
        .def("primes",  &PyOmniFieldPersistence::primes,    "primes over which the matrix differs from the rest")
        .def("column",  [](const PyOmniFieldPersistence& ofp, Index i, BaseElement p)
                        {
                            auto it = ofp.zp_chains().find(i);
                            if (it != ofp.zp_chains().end())
                            {
                                auto pit = it->second.find(p);
                                if (pit != it->second.end())
                                    return pit->second;
                            }
                            return ofp.convert(ofp.q_chains()[i], ofp.zp(p));
                        },                                  "get the column over a specific prime")
        .def("special", &PyOmniFieldPersistence::special,   "test whether the column has a special value over the given prime")
        .def("__len__", &PyOmniFieldPersistence::size,      "size of the persistence object")
        .def("__repr__",    [](const PyOmniFieldPersistence& ofp)
                            { std::ostringstream oss; oss << "OmniFieldPersistence with " << ofp.size() << " columns"; return oss.str(); })
    ;
}
