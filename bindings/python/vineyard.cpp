#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sstream>
namespace py = pybind11;

#include <dionysus/vineyard.h>

#include "field.h"

using PyVineyardMatrix = dionysus::VineyardMatrix<PyZpField>;
using PyVineyard = dionysus::Vineyard<PyZpField>;

void init_vineyard(py::module& m)
{
    using namespace pybind11::literals;
    using Column = std::vector<std::tuple<PyVineyardMatrix::FieldElement, PyVineyardMatrix::Index>>;
    using Columns = std::vector<Column>;
    using Index = PyVineyardMatrix::Index;

    auto make_chain = [](Column entries)
    {
        PyVineyardMatrix::Chain chain;
        chain.reserve(entries.size());
        for (const auto& entry : entries)
            chain.emplace_back(PyVineyardMatrix::Entry { std::get<0>(entry), std::get<1>(entry) });
        return chain;
    };

    auto make_column = [](const PyVineyardMatrix::Chain& chain)
    {
        Column entries;
        entries.reserve(chain.size());
        for (const auto& entry : chain)
            entries.emplace_back(entry.element(), entry.index());
        return entries;
    };

    py::class_<PyVineyardMatrix>(m, "VineyardMatrix", "matrix storage for vineyard updates using stable cell ids")
        .def(py::init<PyZpField, size_t>(), "field"_a = PyZpField(2), "size"_a = 0)
        .def("__len__",      &PyVineyardMatrix::size,      "number of cells")
        .def("field",        &PyVineyardMatrix::field,     "field used for coefficients")
        .def("cell_at",      &PyVineyardMatrix::cell_at,   "cell id at a current filtration position")
        .def("position",     &PyVineyardMatrix::position,  "current filtration position of a cell id")
        .def("low",          &PyVineyardMatrix::low,       "cached low cell id of a column")
        .def("pivot",        &PyVineyardMatrix::pivot,     "column id with the given low cell id")
        .def("refresh_low",  &PyVineyardMatrix::refresh_low, "column"_a,
                                                           "recompute the cached low for one column under the current order")
        .def("set_low",      &PyVineyardMatrix::set_low,   "column"_a, "low"_a,
                                                           "set the cached low for one column")
        .def("contains",     &PyVineyardMatrix::contains,  "column"_a, "row"_a,
                                                           "test whether a stable-id-sorted column contains a row id")
        .def("transpose_position", &PyVineyardMatrix::transpose_position, "position"_a,
                                                            "transpose adjacent filtration positions and return the swapped stable cell ids")
        .def("set_column",   [make_chain](PyVineyardMatrix& vm, Index column, Column entries)
                               {
                                   vm.set_column(column, make_chain(std::move(entries)));
                               }, "column"_a, "entries"_a,
                               "set a column from (coefficient, stable_row_id) entries")
        .def("column",       [make_column](const PyVineyardMatrix& vm, Index column)
                               {
                                   return make_column(vm[column]);
                               }, "column"_a,
                               "return a column as (coefficient, stable_row_id) entries")
        .def("__getitem__",  [make_column](const PyVineyardMatrix& vm, Index column)
                               {
                                   return make_column(vm[column]);
                               }, "column"_a,
                               "return a column as (coefficient, stable_row_id) entries")
        .def_property_readonly("unpaired", [](const PyVineyardMatrix&) { return PyVineyardMatrix::unpaired(); },
                              "index representing lack of pair")
        .def("__repr__",     [](const PyVineyardMatrix& vm)
                              {
                                  std::ostringstream oss;
                                  oss << "VineyardMatrix with " << vm.size() << " cells over " << py::repr(py::cast(vm.field()));
                                  return oss.str();
                               })
    ;

    py::class_<PyVineyard>(m, "Vineyard", "matrix_v-based vineyard state for adjacent transpositions")
        .def(py::init([make_chain](PyZpField field, Columns columns)
                      {
                          PyVineyard::Chains chains;
                          chains.reserve(columns.size());
                          for (auto& column : columns)
                              chains.emplace_back(make_chain(std::move(column)));
                          return new PyVineyard(field, std::move(chains));
                      }), "field"_a = PyZpField(2), "columns"_a = Columns())
        .def("__len__",      &PyVineyard::size,      "number of cells")
        .def("field",        &PyVineyard::field,     "field used for coefficients")
        .def("cell_at",      &PyVineyard::cell_at,   "cell id at a current filtration position")
        .def("position",     &PyVineyard::position,  "current filtration position of a cell id")
        .def("low",          &PyVineyard::low,       "cached low cell id of a reduced column")
        .def("pivot",        &PyVineyard::pivot,     "column id with the given low cell id")
        .def("pair",         &PyVineyard::pair,      "persistence pair of the given cell id")
        .def("transpose_position", &PyVineyard::transpose_position, "position"_a,
                                                            "repair the vineyard state after transposing adjacent filtration positions")
        .def("boundary_column", [make_column](const PyVineyard& v, Index column)
                                {
                                    return make_column(v.boundary_column(column));
                                }, "column"_a,
                                "return an original boundary column as (coefficient, stable_row_id) entries")
        .def("reduced_column", [make_column](const PyVineyard& v, Index column)
                                {
                                    return make_column(v.reduced_column(column));
                                }, "column"_a,
                                "return a reduced column as (coefficient, stable_row_id) entries")
        .def("chain",        [make_column](const PyVineyard& v, Index column)
                                {
                                    return make_column(v.chain(column));
                                }, "column"_a,
                                "return a V-matrix chain as (coefficient, stable_column_id) entries")
        .def_property_readonly("unpaired", [](const PyVineyard&) { return PyVineyard::unpaired(); },
                              "index representing lack of pair")
        .def("__repr__",     [](const PyVineyard& v)
                              {
                                  std::ostringstream oss;
                                  oss << "Vineyard with " << v.size() << " cells over " << py::repr(py::cast(v.field()));
                                  return oss.str();
                              })
    ;
}
