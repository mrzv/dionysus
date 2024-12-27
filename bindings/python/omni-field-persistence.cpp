#include <boost/range/adaptors.hpp>
namespace ba = boost::adaptors;

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "filtration.h"
#include "persistence.h"        // for PyMatrixFiltration
#include "diagram.h"
#include "omni-field-persistence.h"
#include "chain.h"

template<class Filtration>
PyOmniFieldPersistence
omnifield_homology_persistence(const Filtration& filtration)
{
    PyOmniFieldPersistence persistence;
    size_t i = 0;
    for(auto& s : filtration)
    {
        using CellEntry = typename Filtration::Cell::template Entry<PyOmniFieldPersistence::Field>;
        using ChainEntry = dionysus::ChainEntry<PyOmniFieldPersistence::Field, PyOmniFieldPersistence::Index>;
        persistence.add(s.boundary(persistence.field()) |
                                                 ba::transformed([&filtration,i](const CellEntry& e)
                                                 { return ChainEntry(e.element(), filtration.index(e.index(), i)); }));
        ++i;
    }
    return persistence;
}

template<class Filtration>
std::vector<PyDiagram>
py_init_omni_diagrams(const PyOmniFieldPersistence& persistence, const Filtration& f, PyOmniFieldPersistence::BaseElement p)
{
    using Index = PyOmniFieldPersistence::Index;
    return init_diagrams(prime_adapter(persistence, p), f,
                         [](const typename Filtration::Cell& s, Index)  { return s.data(); },        // value
                         [](Index i) -> PyIndex                         { return i; });              // data
}

PYBIND11_MAKE_OPAQUE(PyOmniFieldPersistence::ZpChain);      // persistence.cpp provides a binding for Chain, which is exactly what this is
PYBIND11_MAKE_OPAQUE(PyOmniFieldPersistence::QChain);       // we want to provide our own binding for QChain

void init_omnifield_persistence(py::module& m)
{
    using namespace pybind11::literals;
    m.def("omnifield_homology_persistence",   &omnifield_homology_persistence<PyFiltration>, "filtration"_a,
          "compute homology persistence of the filtration (pair simplices) over all fields at once");
    m.def("omnifield_homology_persistence",   &omnifield_homology_persistence<PyMatrixFiltration>, "filtration"_a,
          "compute homology persistence of the matrix filtration over all fields at once");

    m.def("init_diagrams",      &py_init_omni_diagrams<PyFiltration>,       "ofp"_a, "f"_a, "p"_a,  "initialize diagrams for a specific prime from omnifield persistence and filtration");
    m.def("init_diagrams",      &py_init_omni_diagrams<PyMatrixFiltration>, "ofp"_a, "f"_a, "p"_a,  "initialize diagrams for a specific prime from omnifield persistence and filtration");

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
        .def("column",  [](const PyOmniFieldPersistence& ofp, Index i)
                        {
                            return ofp.q_chains()[i];
                        },                                  "get the column over rationals")
        .def("specials",[](const PyOmniFieldPersistence& ofp)
                        {
                            return ofp.specials();
                        },                                  "get dictionary of special columns mapping to primes")
        .def("special", [](const PyOmniFieldPersistence& ofp, Index i, BaseElement p)
                        {
                            return ofp.special(i,p);
                        },                                  "test whether the column has a special value over the given prime")
        .def("__len__", &PyOmniFieldPersistence::size,      "size of the persistence object")
        .def("__repr__",    [](const PyOmniFieldPersistence& ofp)
                            { std::ostringstream oss; oss << "OmniFieldPersistence with " << ofp.size() << " columns"; return oss.str(); })
    ;

    init_chain<PyOmniFieldPersistence::QChain>(m, "Q");
}
