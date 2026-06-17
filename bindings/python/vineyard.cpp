#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <cmath>
#include <limits>
#include <map>
#include <queue>
#include <sstream>
#include <string>
#include <tuple>
namespace py = pybind11;

#include <dionysus/boundary-matrix.h>
#include <dionysus/vineyard-linear-homotopy.h>
#include <dionysus/vineyard.h>

#include "field.h"
#include "filtration.h"

using PyVineyardMatrix = dionysus::VineyardMatrix<PyZpField>;
using PyVineyardV = dionysus::VineyardV<PyZpField>;
using PyVineyardU = dionysus::VineyardU<PyZpField>;

namespace
{
using Index = PyVineyardMatrix::Index;
using Chain = PyVineyardMatrix::Chain;
using Chains = PyVineyardMatrix::Chains;

using VineyardLinearHomotopyEvent = dionysus::VineyardLinearHomotopyEvent<Index>;
using VineyardSegment = dionysus::VineyardLinearSegment<Index>;
using VineyardVine = dionysus::VineyardLinearVine<Index>;

struct VineyardLinearHomotopyResult
{
    py::object                         vineyard;
    std::vector<VineyardLinearHomotopyEvent> events;
    std::vector<VineyardVine>          vines;
    std::vector<Index>                 final_order;
};

using LinearHomotopyData = dionysus::VineyardLinearHomotopyData<Chains>;

using CrossingCandidate = dionysus::VineyardLinearCrossingCandidate<Index>;
using EventQueue = std::priority_queue<CrossingCandidate>;
using Feature = std::tuple<Index, Index>;

using ActiveVine = dionysus::VineyardLinearActiveVine<Index>;

using ActiveVines = std::map<Feature, ActiveVine>;

template<class Vineyard>
VineyardLinearHomotopyResult run_linear_homotopy(const PyFiltration& filtration,
                                                 const std::vector<double>& values0,
                                                 const std::vector<double>& values1,
                                                 const PyZpField& field)
{
    LinearHomotopyData data = dionysus::prepare_vineyard_linear_homotopy_data<Chains>(filtration, values0, values1, field);
    auto* vineyard = new Vineyard(field, data.boundary);

    VineyardLinearHomotopyResult result;
    ActiveVines active_vines;
    for (const auto& feature : dionysus::vineyard_linear_features<Vineyard, Feature>(*vineyard))
        dionysus::open_vineyard_linear_feature(result, active_vines, feature, 0.0, -1);

    double current_t = 0.0;
    EventQueue event_queue = dionysus::build_vineyard_linear_event_queue<EventQueue>(*vineyard, data, current_t);

    while (current_t < 1.0)
    {
        CrossingCandidate candidate;
        Index position = 0;
        if (!dionysus::pop_next_vineyard_linear_candidate(event_queue, *vineyard, data, current_t, candidate, position) ||
            candidate.time >= 1.0)
        {
            current_t = 1.0;
            break;
        }

        current_t = candidate.time;
        Index first = vineyard->cell_at(position);
        Index second = vineyard->cell_at(position + 1);
        dionysus::validate_vineyard_linear_transposition(data.boundary, first, second);
        dionysus::record_vineyard_linear_transposition<VineyardLinearHomotopyResult, ActiveVines, Vineyard, LinearHomotopyData, Feature>(result, active_vines, *vineyard, data, current_t, position);
        dionysus::push_vineyard_linear_candidate_neighborhood(event_queue, *vineyard, data, current_t, position);
    }

    std::vector<Index> expected_final_order;
    expected_final_order.reserve(filtration.size());
    for (auto input : dionysus::vineyard_linear_endpoint_order(filtration, values1))
        expected_final_order.push_back(data.input_to_stable[input]);

    dionysus::complete_vineyard_linear_endpoint_order<VineyardLinearHomotopyResult, ActiveVines, Vineyard, LinearHomotopyData, Feature>(result, active_vines, *vineyard, data, expected_final_order, 1.0);
    dionysus::close_all_vineyard_linear_features(result, active_vines, data.values0, data.values1, 1.0, -1, Vineyard::unpaired());
    result.final_order = dionysus::current_vineyard_linear_order(*vineyard);
    result.vineyard = py::cast(vineyard, py::return_value_policy::take_ownership);
    return result;
}

}       // namespace

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

    auto make_chains = [make_chain](Columns columns)
    {
        PyVineyardMatrix::Chains chains;
        chains.reserve(columns.size());
        for (auto& column : columns)
            chains.emplace_back(make_chain(std::move(column)));
        return chains;
    };

    py::class_<VineyardLinearHomotopyEvent>(m, "VineyardLinearHomotopyEvent", "one adjacent transposition in a vineyard linear homotopy")
        .def_readonly("time",               &VineyardLinearHomotopyEvent::time)
        .def_readonly("position",           &VineyardLinearHomotopyEvent::position)
        .def_readonly("first",              &VineyardLinearHomotopyEvent::first)
        .def_readonly("second",             &VineyardLinearHomotopyEvent::second)
        .def_readonly("first_pair_before",  &VineyardLinearHomotopyEvent::first_pair_before)
        .def_readonly("second_pair_before", &VineyardLinearHomotopyEvent::second_pair_before)
        .def_readonly("first_pair_after",   &VineyardLinearHomotopyEvent::first_pair_after)
        .def_readonly("second_pair_after",  &VineyardLinearHomotopyEvent::second_pair_after)
        .def_readonly("pairing_switched",   &VineyardLinearHomotopyEvent::pairing_switched)
    ;

    py::class_<VineyardSegment>(m, "VineyardSegment", "one linear segment of a persistence vine")
        .def_readonly("t0",         &VineyardSegment::t0)
        .def_readonly("t1",         &VineyardSegment::t1)
        .def_readonly("birth0",     &VineyardSegment::birth0)
        .def_readonly("birth1",     &VineyardSegment::birth1)
        .def_readonly("death0",     &VineyardSegment::death0)
        .def_readonly("death1",     &VineyardSegment::death1)
        .def_readonly("birth_cell", &VineyardSegment::birth_cell)
        .def_readonly("death_cell", &VineyardSegment::death_cell)
        .def_readonly("event0",     &VineyardSegment::event0)
        .def_readonly("event1",     &VineyardSegment::event1)
    ;

    py::class_<VineyardVine>(m, "VineyardVine", "piecewise-linear persistence vine")
        .def_readonly("segments", &VineyardVine::segments)
    ;

    py::class_<VineyardLinearHomotopyResult>(m, "VineyardLinearHomotopyResult", "result of a vineyard linear homotopy")
        .def_readonly("vineyard",    &VineyardLinearHomotopyResult::vineyard)
        .def_readonly("events",      &VineyardLinearHomotopyResult::events)
        .def_readonly("vines",       &VineyardLinearHomotopyResult::vines)
        .def_readonly("final_order", &VineyardLinearHomotopyResult::final_order)
    ;

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
        .def("transpose", &PyVineyardMatrix::transpose, "position"_a,
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

    py::class_<PyVineyardV>(m, "VineyardV", "matrix_v-based vineyard state for adjacent transpositions")
        .def(py::init([make_chains](Columns columns, PyZpField field)
                        {
                            return new PyVineyardV(field, make_chains(std::move(columns)));
                        }), "columns"_a = Columns(), "field"_a = PyZpField(2))
        .def(py::init([](const PyFiltration& filtration, PyZpField field)
                       {
                           return new PyVineyardV(field, dionysus::make_boundary_chains<Chains>(filtration, field));
                       }), "filtration"_a, "field"_a = PyZpField(2))
        .def("__len__",      &PyVineyardV::size,      "number of cells")
        .def("field",        &PyVineyardV::field,     "field used for coefficients")
        .def("cell_at",      &PyVineyardV::cell_at,   "cell id at a current filtration position")
        .def("position",     &PyVineyardV::position,  "current filtration position of a cell id")
        .def("low",          &PyVineyardV::low,       "cached low cell id of a reduced column")
        .def("pivot",        &PyVineyardV::pivot,     "column id with the given low cell id")
        .def("pair",         &PyVineyardV::pair,      "persistence pair of the given cell id")
        .def("transpose", &PyVineyardV::transpose, "position"_a,
                                                               "repair the vineyard state after transposing adjacent filtration positions and return the swapped stable cell ids plus whether the pairing switched")
        .def("reduced_column", [make_column](const PyVineyardV& v, Index column)
                                 {
                                     return make_column(v.reduced_column(column));
                                 }, "column"_a,
                                 "return a reduced column as (coefficient, stable_row_id) entries")
        .def("chain",        [make_column](const PyVineyardV& v, Index column)
                                 {
                                     return make_column(v.basis(column));
                                 }, "column"_a,
                                 "return a V-matrix chain as (coefficient, stable_column_id) entries")
        .def_property_readonly("unpaired", [](const PyVineyardV&) { return PyVineyardV::unpaired(); },
                               "index representing lack of pair")
        .def("__repr__",     [](const PyVineyardV& v)
                               {
                                   std::ostringstream oss;
                                   oss << "VineyardV with " << v.size() << " cells over " << py::repr(py::cast(v.field()));
                                   return oss.str();
                               })
    ;

    py::class_<PyVineyardU>(m, "VineyardU", "matrix_u-based vineyard state for adjacent transpositions")
        .def(py::init([make_chains](Columns columns, PyZpField field)
                        {
                            return new PyVineyardU(field, make_chains(std::move(columns)));
                        }), "columns"_a = Columns(), "field"_a = PyZpField(2))
        .def(py::init([](const PyFiltration& filtration, PyZpField field)
                       {
                           return new PyVineyardU(field, dionysus::make_boundary_chains<Chains>(filtration, field));
                       }), "filtration"_a, "field"_a = PyZpField(2))
        .def("__len__",      &PyVineyardU::size,      "number of cells")
        .def("field",        &PyVineyardU::field,     "field used for coefficients")
        .def("cell_at",      &PyVineyardU::cell_at,   "cell id at a current filtration position")
        .def("position",     &PyVineyardU::position,  "current filtration position of a cell id")
        .def("low",          &PyVineyardU::low,       "cached low cell id of a reduced column")
        .def("pivot",        &PyVineyardU::pivot,     "column id with the given low cell id")
        .def("pair",         &PyVineyardU::pair,      "persistence pair of the given cell id")
        .def("transpose", &PyVineyardU::transpose, "position"_a,
                                                               "repair the vineyard state after transposing adjacent filtration positions and return the swapped stable cell ids plus whether the pairing switched")
        .def("reduced_column", [make_column](const PyVineyardU& v, Index column)
                                 {
                                     return make_column(v.reduced_column(column));
                                 }, "column"_a,
                                 "return a reduced column as (coefficient, stable_row_id) entries")
        .def("trail",        [make_column](const PyVineyardU& v, Index row)
                                 {
                                     return make_column(v.basis(row));
                                 }, "row"_a,
                                 "return a U-matrix trail row as (coefficient, stable_column_id) entries")
        .def_property_readonly("unpaired", [](const PyVineyardU&) { return PyVineyardU::unpaired(); },
                               "index representing lack of pair")
        .def("__repr__",     [](const PyVineyardU& v)
                               {
                                   std::ostringstream oss;
                                   oss << "VineyardU with " << v.size() << " cells over " << py::repr(py::cast(v.field()));
                                   return oss.str();
                               })
    ;

    m.def("Vineyard", [make_chains](Columns columns, PyZpField field, std::string method) -> py::object
          {
              if (method == "matrix_v")
                  return py::cast(new PyVineyardV(field, make_chains(std::move(columns))), py::return_value_policy::take_ownership);
              if (method == "matrix_u")
                  return py::cast(new PyVineyardU(field, make_chains(std::move(columns))), py::return_value_policy::take_ownership);
              throw py::value_error("unknown vineyard method: " + method);
          }, "columns"_a = Columns(), "field"_a = PyZpField(2), "method"_a = "matrix_v",
          "construct a vineyard state using method='matrix_v' or method='matrix_u'");

    m.def("Vineyard", [](const PyFiltration& filtration, PyZpField field, std::string method) -> py::object
          {
              if (method == "matrix_v")
                  return py::cast(new PyVineyardV(field, dionysus::make_boundary_chains<Chains>(filtration, field)), py::return_value_policy::take_ownership);
              if (method == "matrix_u")
                  return py::cast(new PyVineyardU(field, dionysus::make_boundary_chains<Chains>(filtration, field)), py::return_value_policy::take_ownership);
              throw py::value_error("unknown vineyard method: " + method);
          }, "filtration"_a, "field"_a = PyZpField(2), "method"_a = "matrix_v",
          "construct a vineyard state using method='matrix_v' or method='matrix_u'");

    m.def("vineyard_linear_homotopy", [](const PyFiltration& filtration,
                                           std::vector<double> values0,
                                           std::vector<double> values1,
                                           PyZpField field,
                                           std::string method)
          {
              if (method == "matrix_v")
                  return run_linear_homotopy<PyVineyardV>(filtration, values0, values1, field);
              if (method == "matrix_u")
                  return run_linear_homotopy<PyVineyardU>(filtration, values0, values1, field);
              throw py::value_error("unknown vineyard method: " + method);
           }, "filtration"_a, "values0"_a, "values1"_a, "field"_a = PyZpField(2), "method"_a = "matrix_v",
           "compute a vineyard linear homotopy between two strict filtration functions while preserving vine identity combinatorially");
}
