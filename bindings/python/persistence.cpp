#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <boost/range/adaptors.hpp>
namespace ba = boost::adaptors;

#include <dionysus/row-reduction.h>
#include <dionysus/ordinary-persistence.h>
#include <dionysus/standard-reduction.h>
#include <dionysus/clearing-reduction.h>

#include "field.h"
#include "filtration.h"
#include "persistence.h"
#include "diagram.h"
#include "chain.h"
#include "progress.h"

template<class Filtration, class Relative>
py::object
compute_homology_persistence(const Filtration& filtration, const Relative& relative, PyZpField::Element prime, std::string method, const Progress& progress)
{
    PyZpField field(prime);

    if (method == "clearing")
    {
        using Persistence = dionysus::OrdinaryPersistence<PyZpField>;
        using Reduction   = dionysus::ClearingReduction<Persistence>;
        Persistence persistence(field);
        Reduction   reduce(persistence);
        reduce(filtration, relative, &Reduction::no_report_pair, progress);
        return py::cast(std::move(reduce.persistence()));
    }
    else if (method == "row")
    {
        using Reduction = dionysus::RowReduction<PyZpField>;
        Reduction   reduce(field);
        reduce(filtration, relative, &Reduction::no_report_pair, progress);
        return py::cast(std::move(reduce.persistence()));
    } else if (method == "column")
    {
        using Persistence = dionysus::OrdinaryPersistence<PyZpField>;
        using Reduction   = dionysus::StandardReduction<Persistence>;
        Persistence persistence(field);
        Reduction   reduce(persistence);
        reduce(filtration, relative, &Reduction::no_report_pair, progress);
        return py::cast(std::move(reduce.persistence()));
    } else if (method == "column_no_negative")
    {
        using Persistence = dionysus::OrdinaryPersistenceNoNegative<PyZpField>;
        using Reduction   = dionysus::StandardReduction<Persistence>;
        Persistence persistence(field);
        Reduction   reduce(persistence);
        reduce(filtration, relative, &Reduction::no_report_pair, progress);
        return py::cast(std::move(reduce.persistence()));
    } else if (method == "matrix_v")
    {
        using Persistence = dionysus::OrdinaryPersistenceWithV<PyZpField>;
        using Reduction   = dionysus::StandardReduction<Persistence>;
        Persistence persistence(field);
        Reduction   reduce(persistence);
        reduce(filtration, relative, &Reduction::no_report_pair, progress);
        auto v = std::move(reduce.persistence().visitor<0>().v_);
        auto r = std::move(reduce.persistence());
        return py::cast(std::make_tuple(r, v));
    } else if (method == "matrix_v_no_negative")
    {
        using Persistence = dionysus::OrdinaryPersistenceNoNegativeWithV<PyZpField>;
        using Reduction   = dionysus::StandardReduction<Persistence>;
        Persistence persistence(field);
        Reduction   reduce(persistence);
        reduce(filtration, relative, &Reduction::no_report_pair, progress);
        auto v = std::move(reduce.persistence().visitor<1>().v_);
        auto r = std::move(reduce.persistence());
        return py::cast(std::make_tuple(r,v));
    } else
        throw std::runtime_error("Unknown method: " + method);
}

template<class Filtration>
py::object
homology_persistence(const Filtration& filtration, PyZpField::Element prime, std::string method, bool progress)
{
    using Cell = typename Filtration::Cell;
    if (progress)
        return compute_homology_persistence(filtration, [](const Cell&) { return false; }, prime, method, ShowProgress(filtration.size()));
    else
        return compute_homology_persistence(filtration, [](const Cell&) { return false; }, prime, method, NoProgress());
}

py::object
relative_homology_persistence(const PyFiltration& filtration, const PyFiltration& relative, PyZpField::Element prime, std::string method, bool progress)
{
    using Cell = PyFiltration::Cell;
    if (progress)
        return compute_homology_persistence(filtration, [&relative](const Cell& c) { return relative.contains(c); }, prime, method, ShowProgress(filtration.size()));
    else
        return compute_homology_persistence(filtration, [&relative](const Cell& c) { return relative.contains(c); }, prime, method, NoProgress());
}

template<class PyReducedMatrix, class Filtration>
std::vector<PyDiagram>
py_init_diagrams(const PyReducedMatrix& m, const Filtration& f)
{
    using Index = typename PyReducedMatrix::Index;
    return init_diagrams(m, f,
                         [](const typename Filtration::Cell& s, Index)      { return s.data(); },       // value
                         [](Index i) -> PyIndex                             { return i; });             // data
}

bool
homologous(PyReducedMatrix& m, PyReducedMatrix::Chain z1, PyReducedMatrix::Chain z2)
{
    using Entry = PyReducedMatrix::Entry;
    auto entry_cmp = [&m](const Entry& e1, const Entry& e2) { return m.cmp()(e1.index(), e2.index()); };

    std::sort(z1.begin(), z1.end(), entry_cmp);
    std::sort(z2.begin(), z2.end(), entry_cmp);

    // z1 -= z2
    dionysus::Chain<PyReducedMatrix::Chain>::addto(z1, m.field().neg(m.field().id()), z2, m.field(), entry_cmp);
    m.reduce(m.size(), z1);        // dummy m.size() index (only used for callbacks, which are empty here, so unused)
    return z1.empty();
}

PYBIND11_MAKE_OPAQUE(PyReducedMatrix::Chain);      // we want to provide our own binding for Chain
PYBIND11_MAKE_OPAQUE(PyMatrixFiltration::Cell::BoundaryChain<>);      // we want to provide our own binding for BoundaryChain

template<class PyReducedMatrix>
void export_reduced_matrix(py::module& m, std::string name)
{
    using namespace pybind11::literals;

    // for pickling
    using Column = std::vector<std::tuple<typename PyReducedMatrix::FieldElement, typename PyReducedMatrix::Index>>;
    using Columns = std::vector<Column>;
    using Pairs = std::vector<typename PyReducedMatrix::Index>;
    using Skips = std::vector<bool>;
    using Index = typename PyReducedMatrix::Index;
    using Chain = typename PyReducedMatrix::Chain;

    py::class_<PyReducedMatrix>(m, name.c_str(), "matrix, where each column has a lowest non-zero entry in a unique row; supports iteration and indexing")
        .def(py::init([](PyZpField field, size_t sz)
                      {
                        PyReducedMatrix* m = new PyReducedMatrix(field);
                        m->resize(sz);
                        return m;
                      }), "field"_a = PyZpField(2), "size"_a = 0)
        .def("__len__",     &PyReducedMatrix::size,         "size of the matrix")
        .def("__getitem__", &PyReducedMatrix::operator[],   "access the column at a given index")
        .def("__setitem__", [](PyReducedMatrix* m, Index i, Column c) { m->set(i,std::move(c)); },
                                                            "set the column at a given index")
        .def("__setitem__", [](PyReducedMatrix* m, Index i, Chain c) { m->set(i,std::move(c)); },
                                                            "set the column at a given index")
        .def("pair",        &PyReducedMatrix::pair,         "pair of the given index")
        .def_property_readonly("unpaired",      [](const PyReducedMatrix&) { return PyReducedMatrix::unpaired(); },
                               "index representing lack of pair")
        .def("homologous",  &homologous,                    "test if two cycles are homologous")
        .def("resize",      &PyReducedMatrix::resize,       "resize the number of columns")
        .def("field",       &PyReducedMatrix::field,        "access the field used for the reduction")
        .def("__iter__",    [](const PyReducedMatrix& rm)   { return py::make_iterator(rm.columns().begin(), rm.columns().end()); },
                                py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */,
                                "iterate over the columns of the matrix")
        .def("__repr__",    [](const PyReducedMatrix& rm)
                            { std::ostringstream oss; oss << "Reduced matrix with " << rm.size() << " columns over " << py::repr(py::cast(rm.field())); return oss.str(); })
        .def(py::pickle(
            [](const PyReducedMatrix& m)        // __getstate__
            {
                Columns columns; columns.reserve(m.size());
                Pairs pairs; pairs.reserve(m.size());
                Skips skips; skips.reserve(m.size());

                typename PyReducedMatrix::Index i = 0;
                while (i < m.size())
                {
                    Column c;
                    for (auto& e : m[i])
                        c.emplace_back(e.element(), e.index());
                    columns.emplace_back(std::move(c));
                    pairs.emplace_back(m.pair(i));
                    skips.emplace_back(m.skip(i));
                    ++i;
                }

                return py::make_tuple(m.field().prime(), columns, pairs, skips);
            },
            [](py::tuple t)                     // __setstate__
            {
                if (t.size() != 1)
                    throw std::runtime_error("Invalid state!");

                auto prime = t[0].cast<typename PyReducedMatrix::FieldElement>();
                Columns columns = t[1].cast<Columns>();
                Pairs pairs = t[2].cast<Pairs>();
                Skips skips = t[3].cast<Skips>();
                PyReducedMatrix m(prime);
                m.resize(columns.size());
                typename PyReducedMatrix::Index i = 0;
                for (auto& c : columns)
                {
                    m.set(i, c | ba::transformed([](const typename Column::value_type& e)
                                                   {
                                                     return typename PyReducedMatrix::Entry { std::get<0>(e), std::get<1>(e) };
                                                   }));
                    m.set_pair(i,pairs[i]);
                    m.set_skip(i,skips[i]);
                    ++i;
                }

                return m;
            }
        ));
    ;

    m.def("init_diagrams",      &py_init_diagrams<PyReducedMatrix, PyFiltration>,            "m"_a, "f"_a,  "initialize diagrams from reduced matrix and filtration");
    m.def("init_diagrams",      &py_init_diagrams<PyReducedMatrix, PyMatrixFiltration>,      "m"_a, "f"_a,  "initialize diagrams from reduced matrix and filtration");
    m.def("init_diagrams",      &py_init_diagrams<PyReducedMatrix, PyMultiFiltration>,       "m"_a, "f"_a,  "initialize diagrams from reduced matrix and filtration");
    m.def("init_diagrams",      &py_init_diagrams<PyReducedMatrix, PyLinkedMultiFiltration>, "m"_a, "f"_a,  "initialize diagrams from reduced matrix and filtration");
}

void init_persistence(py::module& m)
{
    using namespace pybind11::literals;
    m.def("homology_persistence",   &homology_persistence<PyFiltration>,
          "filtration"_a, "prime"_a = 2, "method"_a = "clearing", "progress"_a = false,
          "compute homology persistence of the filtration (pair simplices); method is one of `clearing`, `row`, `column`, or `column_no_negative`");
    m.def("homology_persistence",   &homology_persistence<PyMatrixFiltration>,
          "filtration"_a, "prime"_a = 2, "method"_a = "clearing", "progress"_a = false,
          "compute homology persistence of the filtration (pair simplices); method is one of `clearing`, `row`, `column`, or `column_no_negative`");
    m.def("homology_persistence",   &homology_persistence<PyMultiFiltration>,
          "filtration"_a, "prime"_a = 2, "method"_a = "clearing", "progress"_a = false,
          "compute homology persistence of the filtration (pair simplices); method is one of `clearing`, `row`, `column`, or `column_no_negative`");
    m.def("homology_persistence",   &homology_persistence<PyLinkedMultiFiltration>,
          "filtration"_a, "prime"_a = 2, "method"_a = "clearing", "progress"_a = false,
          "compute homology persistence of the filtration (pair simplices); method is one of `clearing`, `row`, `column`, or `column_no_negative`");
    m.def("homology_persistence",   &relative_homology_persistence,
          "filtration"_a, "relative"_a, "prime"_a = 2, "method"_a = "clearing", "progress"_a = false,
          "compute homology persistence of the filtration, relative to a subcomplex; method is one of `clearing`, `row`, `column`, or `column_no_negative`");

    py::class_<PyMatrixFiltration::Cell>(m, "MatrixFiltrationCell", "Cell-like adapter for a matrix column")
        .def("__repr__",    [](const PyMatrixFiltration::Cell& mfc)
                            { std::ostringstream oss; oss << "Cell " << mfc.i(); return oss.str(); })
        .def("dimension",   &PyMatrixFiltration::Cell::dimension, "cell dimension")
        .def("boundary",    [](const PyMatrixFiltration::Cell& mfc) { return mfc.boundary(); },
                            "boundary of the cell (the column in the matrix)")
    ;

    export_reduced_matrix<PyReducedMatrix>(m, "ReducedMatrix");
    export_reduced_matrix<PyReducedMatrixNoNegative>(m, "ReducedMatrixNoNegative");
    export_reduced_matrix<PyReducedMatrixWithV>(m, "ReducedMatrixWithV");
    export_reduced_matrix<PyReducedMatrixNoNegativeWithV>(m, "ReducedMatrixNoNegativeWithV");
    init_chain<typename PyReducedMatrix::Chain>(m);

    py::class_<PyMatrixFiltration>(m, "MatrixFiltration", "adapter to turn ReducedMatrix into something that looks and acts like a filtration")
        .def(py::init<PyReducedMatrix, Dimensions, Values>())
        .def("dimensions",  &PyMatrixFiltration::dimensions,    "list of cell dimensions")
        .def("values",      &PyMatrixFiltration::values,        "list of cell values")
        .def("__len__",     &PyMatrixFiltration::size,          "size of the matrix")
        .def("__getitem__", &PyMatrixFiltration::operator[],    "access the 'cell' (column) at a given index")
        .def("__repr__",    [](const PyMatrixFiltration& mf)
                            { std::ostringstream oss; oss << "MatrixFiltration with " << mf.size() << " cells"; return oss.str(); })
    ;
    init_chain<PyMatrixFiltration::Cell::BoundaryChain<>>(m, "MF");
}
