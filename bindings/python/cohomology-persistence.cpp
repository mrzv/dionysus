#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <dionysus/cohomology-persistence.h>
#include <dionysus/standard-reduction.h>

#include "cohomology-persistence.h"
#include "filtration.h"
#include "chain.h"
#include "diagram.h"

PyCohomologyPersistence
cohomology_persistence(const PyFiltration& filtration, PyZpField::Element prime, bool keep_cocycles)
{
    PyZpField field(prime);

    using Persistence = PyCohomologyPersistence;
    using Reduction   = dionysus::StandardReduction<Persistence>;

    Persistence persistence(field);
    persistence.keep_cocycles = keep_cocycles;

    Reduction   reduce(persistence);
    reduce(filtration);

    return persistence;
}

std::vector<PyDiagram>
py_init_diagrams(const PyCohomologyPersistence& m, const PyFiltration& f)
{
    return init_diagrams(m, f,
                         [](const PySimplex& s)                             { return s.data(); },       // value
                         [](PyCohomologyPersistence::Index i) -> PyIndex    { return i; });             // data
}

PYBIND11_MAKE_OPAQUE(PyCohomologyPersistence::Column);      // we want to provide our own binding for the cochain

void init_cohomology_persistence(py::module& m)
{
    using namespace pybind11::literals;
    m.def("cohomology_persistence",   &cohomology_persistence, "filtration"_a, "prime"_a = 2, "keep_cocycles"_a = false,
          "compute cohomology persistence of the filtration");

    m.def("init_diagrams",      &py_init_diagrams,  "m"_a, "f"_a,  "initialize diagrams from cohomology persistence and filtration");

    py::class_<PyCohomologyPersistence>(m, "CohomologyPersistence", "representation of pairs and alive cocycles")
        .def(py::init<PyZpField>())
        .def("__len__",     &PyCohomologyPersistence::size,         "size of the matrix")
        .def("pair",        &PyCohomologyPersistence::pair,         "pair of the given index")
        .def_property_readonly("unpaired",      [](const PyCohomologyPersistence&) { return PyCohomologyPersistence::unpaired(); },
                               "index representing lack of pair")
        .def("__iter__",    [](const PyCohomologyPersistence& rm)   { return py::make_iterator(rm.columns().begin(), rm.columns().end()); },
                                py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */,
                                "iterate over the column heads of the matrix")
        .def("cocycle",     [](const PyCohomologyPersistence& rm, PyCohomologyPersistence::Index idx)
                            {
                                if (rm.pair(idx) != rm.unpaired())
                                    return rm.chain(idx);
                                else        // search through the alive cocycles
                                    for (auto& x : rm.columns())
                                        if (x.index() == idx)
                                            return x.chain;
                                throw py::index_error();
                            },                                      "cocycle that died with the addition of the cell at the given index")
        .def("__repr__",    [](const PyCohomologyPersistence& rm)
                            { std::ostringstream oss; oss << "Cohomology persistence of " << rm.size() << " cells"; return oss.str(); })
    ;

    using ColumnHead = PyCohomologyPersistence::ColumnHead;
    py::class_<ColumnHead>(m, "CohomologyPersistenceColumnHead")
        .def_property_readonly("index",     [](const ColumnHead& x) { return x.index(); },
                                            "index when the column was added")
        .def_readonly("cocycle",            &ColumnHead::chain,   "still alive cocycle")
        .def("__repr__",                    [](const ColumnHead& x)
                                            { std::ostringstream oss; oss << "cocycle added at index " << x.index(); return oss.str(); })
    ;

    init_chain<PyCohomologyPersistence::Column>(m, "Co");
}
