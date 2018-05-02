#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
namespace py = pybind11;

#include <dionysus/row-reduction.h>
#include <dionysus/ordinary-persistence.h>
#include <dionysus/standard-reduction.h>
#include <dionysus/clearing-reduction.h>

#include <dionysus/dlog/progress.h>

#include "field.h"
#include "filtration.h"
#include "persistence.h"
#include "diagram.h"
#include "chain.h"

struct ShowProgress
{
        ShowProgress(size_t total):
            progress(total)                     {}

    void    operator()() const                  { ++progress; }

    mutable dlog::progress  progress;
};

struct NoProgress
{
    void    operator()() const                  {}
};

template<class Relative, class Progress>
PyReducedMatrix
compute_homology_persistence(const PyFiltration& filtration, const Relative& relative, PyZpField::Element prime, std::string method, const Progress& progress)
{
    PyZpField field(prime);

    py::scoped_ostream_redirect stream(
        std::cout,                               // std::ostream&
        py::module::import("sys").attr("stdout") // Python output
    );

    if (method == "clearing")
    {
        using Persistence = dionysus::OrdinaryPersistence<PyZpField>;
        using Reduction   = dionysus::ClearingReduction<Persistence>;
        Persistence persistence(field);
        Reduction   reduce(persistence);
        reduce(filtration, relative, &Reduction::no_report_pair, progress);
        return std::move(reduce.persistence());
    }
    else if (method == "row")
    {
        using Reduction = dionysus::RowReduction<PyZpField>;
        Reduction   reduce(field);
        reduce(filtration, relative, &Reduction::no_report_pair, progress);
        return reduce.persistence();
    } else if (method == "column")
    {
        using Persistence = dionysus::OrdinaryPersistence<PyZpField>;
        using Reduction   = dionysus::StandardReduction<Persistence>;
        Persistence persistence(field);
        Reduction   reduce(persistence);
        reduce(filtration, relative, &Reduction::no_report_pair, progress);
        return std::move(reduce.persistence());
    } else if (method == "column_no_negative")
    {
        using Persistence = dionysus::OrdinaryPersistenceNoNegative<PyZpField>;
        using Reduction   = dionysus::StandardReduction<Persistence>;
        Persistence persistence(field);
        Reduction   reduce(persistence);
        reduce(filtration, relative, &Reduction::no_report_pair, progress);
        return std::move(reduce.persistence());
    } else
        throw std::runtime_error("Unknown method: " + method);
}

PyReducedMatrix
homology_persistence(const PyFiltration& filtration, PyZpField::Element prime, std::string method, bool progress)
{
    using Cell = PyFiltration::Cell;
    if (progress)
        return compute_homology_persistence(filtration, [](const Cell&) { return false; }, prime, method, ShowProgress(filtration.size()));
    else
        return compute_homology_persistence(filtration, [](const Cell&) { return false; }, prime, method, NoProgress());
}

PyReducedMatrix
relative_homology_persistence(const PyFiltration& filtration, const PyFiltration& relative, PyZpField::Element prime, std::string method, bool progress)
{
    using Cell = PyFiltration::Cell;
    if (progress)
        return compute_homology_persistence(filtration, [&relative](const Cell& c) { return relative.contains(c); }, prime, method, ShowProgress(filtration.size()));
    else
        return compute_homology_persistence(filtration, [&relative](const Cell& c) { return relative.contains(c); }, prime, method, NoProgress());
}

std::vector<PyDiagram>
py_init_diagrams(const PyReducedMatrix& m, const PyFiltration& f)
{
    return init_diagrams(m, f,
                         [](const PySimplex& s)                     { return s.data(); },       // value
                         [](PyReducedMatrix::Index i) -> PyIndex    { return i; });             // data
}

PYBIND11_MAKE_OPAQUE(PyReducedMatrix::Chain);      // we want to provide our own binding for Chain

void init_persistence(py::module& m)
{
    using namespace pybind11::literals;
    m.def("homology_persistence",   &homology_persistence,
          "filtration"_a, "prime"_a = 2, "method"_a = "clearing", "progress"_a = false,
          "compute homology persistence of the filtration (pair simplices); method is one of `clearing`, `row`, `column`, or `column_no_negative`");
    m.def("homology_persistence",   &relative_homology_persistence,
          "filtration"_a, "relative"_a, "prime"_a = 2, "method"_a = "clearing", "progress"_a = false,
          "compute homology persistence of the filtration, relative to a subcomplex; method is one of `clearing`, `row`, `column`, or `column_no_negative`");

    m.def("init_diagrams",      &py_init_diagrams,  "m"_a, "f"_a,  "initialize diagrams from reduced matrix and filtration");

    py::class_<PyReducedMatrix>(m, "ReducedMatrix", "matrix, where each column has a lowest non-zero entry in a unique row; supports iteration and indexing")
        .def(py::init<PyZpField>())
        .def("__len__",     &PyReducedMatrix::size,         "size of the matrix")
        .def("__getitem__", &PyReducedMatrix::operator[],   "access the column at a given index")
        .def("pair",        &PyReducedMatrix::pair,         "pair of the given index")
        .def_property_readonly("unpaired",      [](const PyReducedMatrix&) { return PyReducedMatrix::unpaired(); },
                               "index representing lack of pair")
        .def("__iter__",    [](const PyReducedMatrix& rm)   { return py::make_iterator(rm.columns().begin(), rm.columns().end()); },
                                py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */,
                                "iterate over the columns of the matrix")
        .def("__repr__",    [](const PyReducedMatrix& rm)
                            { std::ostringstream oss; oss << "Reduced matrix with " << rm.size() << " columns"; return oss.str(); })
    ;

    init_chain<PyReducedMatrix::Chain>(m);
}
