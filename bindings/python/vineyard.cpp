#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sstream>
namespace py = pybind11;

#include <dionysus/vineyard.h>

#include "field.h"

using PyVineyardMatrix = dionysus::VineyardMatrix<PyZpField>;

void init_vineyard(py::module& m)
{
    using namespace pybind11::literals;
    using Column = std::vector<std::tuple<PyVineyardMatrix::FieldElement, PyVineyardMatrix::Index>>;
    using Index = PyVineyardMatrix::Index;

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
        .def("set_column",   [](PyVineyardMatrix& vm, Index column, Column entries)
                              {
                                  PyVineyardMatrix::Chain chain;
                                  chain.reserve(entries.size());
                                  for (const auto& entry : entries)
                                      chain.emplace_back(PyVineyardMatrix::Entry { std::get<0>(entry), std::get<1>(entry) });
                                  vm.set_column(column, std::move(chain));
                              }, "column"_a, "entries"_a,
                              "set a column from (coefficient, stable_row_id) entries")
        .def("column",       [](const PyVineyardMatrix& vm, Index column)
                              {
                                  Column entries;
                                  for (const auto& entry : vm[column])
                                      entries.emplace_back(entry.element(), entry.index());
                                  return entries;
                              }, "column"_a,
                              "return a column as (coefficient, stable_row_id) entries")
        .def("__getitem__",  [](const PyVineyardMatrix& vm, Index column)
                              {
                                  Column entries;
                                  for (const auto& entry : vm[column])
                                      entries.emplace_back(entry.element(), entry.index());
                                  return entries;
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
}
