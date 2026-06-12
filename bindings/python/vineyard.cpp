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

constexpr double epsilon = 1e-10;

struct VineyardLinearHomotopyEvent
{
    double  time;
    Index   position;
    Index   first;
    Index   second;
    Index   first_pair_before;
    Index   second_pair_before;
    Index   first_pair_after;
    Index   second_pair_after;
};

struct VineyardSegment
{
    double  t0;
    double  t1;
    double  birth0;
    double  birth1;
    double  death0;
    double  death1;
    Index   birth_cell;
    Index   death_cell;
    int     event0;
    int     event1;
};

struct VineyardVine
{
    std::vector<VineyardSegment> segments;
};

struct VineyardLinearHomotopyResult
{
    py::object                         vineyard;
    std::vector<VineyardLinearHomotopyEvent> events;
    std::vector<VineyardVine>          vines;
    std::vector<Index>                 final_order;
};

struct LinearHomotopyData
{
    Chains                  boundary;
    std::vector<Index>      stable_to_input;
    std::vector<Index>      input_to_stable;
    std::vector<double>     values0;
    std::vector<double>     values1;
    std::vector<unsigned>   dimensions;
};

struct CrossingCandidate
{
    double  time;
    Index   left;
    Index   right;
    double  priority_slope;
    unsigned priority_dimension;
    Index   priority_cell;

    bool operator<(const CrossingCandidate& other) const
    {
        if (std::abs(time - other.time) > epsilon)
            return time > other.time;
        if (std::abs(priority_slope - other.priority_slope) > epsilon)
            return priority_slope > other.priority_slope;
        if (priority_dimension != other.priority_dimension)
            return priority_dimension > other.priority_dimension;
        return priority_cell > other.priority_cell;
    }
};

using EventQueue = std::priority_queue<CrossingCandidate>;
using Feature = std::tuple<Index, Index>;

struct ActiveVine
{
    size_t  vine_index;
    double  t0;
    int     event0;
    Index   birth;
    Index   death;
};

using ActiveVines = std::map<Feature, ActiveVine>;

double value_at(const std::vector<double>& values0, const std::vector<double>& values1, Index cell, double t)
{
    return values0[cell] + t * (values1[cell] - values0[cell]);
}

double slope(const std::vector<double>& values0, const std::vector<double>& values1, Index cell)
{
    return values1[cell] - values0[cell];
}

template<class Vineyard>
std::vector<Index> current_order(const Vineyard& vineyard)
{
    std::vector<Index> order;
    order.reserve(vineyard.size());
    for (Index p = 0; p < vineyard.size(); ++p)
        order.push_back(vineyard.cell_at(p));
    return order;
}

bool filtration_less(const PyFiltration& filtration,
                     const std::vector<double>& values,
                     Index x,
                     Index y)
{
    if (values[x] != values[y])
        return values[x] < values[y];
    return filtration[x] < filtration[y];
}

void sort_chain(Chain& column)
{
    std::sort(column.begin(), column.end(), [](const PyVineyardMatrix::Entry& x, const PyVineyardMatrix::Entry& y)
              { return x.index() < y.index(); });
}

Chains boundary_from_filtration(const PyFiltration& filtration, const PyZpField& field)
{
    Chains boundary(filtration.size());

    for (Index i = 0; i < filtration.size(); ++i)
    {
        Chain column;
        for (auto it = filtration[i].boundary_begin(field); it != filtration[i].boundary_end(field); ++it)
        {
            Index face = filtration.index(it->index(), filtration.size());
            if (face >= i)
                throw py::value_error("filtration boundary face does not precede simplex");
            column.emplace_back(it->element(), face);
        }
        sort_chain(column);
        boundary[i] = std::move(column);
    }

    return boundary;
}

void validate_endpoint_filtration(const PyFiltration& filtration, const std::vector<double>& values)
{
    for (Index i = 0; i < filtration.size(); ++i)
        for (auto it = filtration[i].boundary_begin(); it != filtration[i].boundary_end(); ++it)
        {
            Index face = filtration.index(*it, filtration.size());
            if (values[face] > values[i] + epsilon)
                throw std::invalid_argument("linear homotopy values are not a filtration on the complex");
        }
}

LinearHomotopyData prepare_linear_homotopy_data(const PyFiltration& filtration,
                                                const std::vector<double>& input_values0,
                                                const std::vector<double>& input_values1,
                                                const PyZpField& field)
{
    if (input_values0.size() != filtration.size() || input_values1.size() != filtration.size())
        throw std::invalid_argument("linear homotopy value vectors must match filtration size");

    validate_endpoint_filtration(filtration, input_values0);
    validate_endpoint_filtration(filtration, input_values1);

    LinearHomotopyData data;
    data.stable_to_input.resize(filtration.size());
    data.input_to_stable.resize(filtration.size());
    for (Index i = 0; i < filtration.size(); ++i)
        data.stable_to_input[i] = i;

    std::sort(data.stable_to_input.begin(), data.stable_to_input.end(),
              [&filtration, &input_values0](Index x, Index y)
              { return filtration_less(filtration, input_values0, x, y); });

    for (Index stable = 0; stable < data.stable_to_input.size(); ++stable)
        data.input_to_stable[data.stable_to_input[stable]] = stable;

    data.values0.resize(filtration.size());
    data.values1.resize(filtration.size());
    data.dimensions.resize(filtration.size());
    data.boundary.resize(filtration.size());

    for (Index stable = 0; stable < data.stable_to_input.size(); ++stable)
    {
        Index input = data.stable_to_input[stable];
        const PySimplex& simplex = filtration[input];
        data.values0[stable] = input_values0[input];
        data.values1[stable] = input_values1[input];
        data.dimensions[stable] = simplex.dimension();

        Chain column;
        for (auto it = simplex.boundary_begin(field); it != simplex.boundary_end(field); ++it)
        {
            Index face_input = filtration.index(it->index(), filtration.size());
            column.emplace_back(it->element(), data.input_to_stable[face_input]);
        }
        sort_chain(column);
        data.boundary[stable] = std::move(column);
    }

    return data;
}

template<class Vineyard>
bool feature_for_cell(const Vineyard& vineyard, Index cell, Feature& feature)
{
    if (cell == Vineyard::unpaired())
        return false;

    Index pair = vineyard.pair(cell);
    if (pair == Vineyard::unpaired())
    {
        if (vineyard.low(cell) == Vineyard::unpaired())
        {
            feature = std::make_tuple(cell, Vineyard::unpaired());
            return true;
        }
        return false;
    }

    if (vineyard.position(cell) < vineyard.position(pair))
        feature = std::make_tuple(cell, pair);
    else
        feature = std::make_tuple(pair, cell);
    return true;
}

template<class Vineyard>
std::vector<Feature> features(const Vineyard& vineyard)
{
    std::vector<Feature> result;
    for (Index cell = 0; cell < vineyard.size(); ++cell)
    {
        Feature feature;
        if (feature_for_cell(vineyard, cell, feature))
            result.push_back(feature);
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

template<class Vineyard>
std::vector<Feature> local_features(const Vineyard& vineyard, const std::vector<Index>& cells)
{
    std::vector<Feature> result;
    for (Index cell : cells)
    {
        Feature feature;
        if (feature_for_cell(vineyard, cell, feature))
            result.push_back(feature);
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

void add_segment(VineyardLinearHomotopyResult& result,
                 const ActiveVine& active,
                 const std::vector<double>& values0,
                 const std::vector<double>& values1,
                 double t1,
                 int event1)
{
    if (t1 <= active.t0 + epsilon)
        return;

    double inf = std::numeric_limits<double>::infinity();
    result.vines[active.vine_index].segments.push_back(VineyardSegment {
        active.t0,
        t1,
        value_at(values0, values1, active.birth, active.t0),
        value_at(values0, values1, active.birth, t1),
        active.death == PyVineyardMatrix::unpaired() ? inf : value_at(values0, values1, active.death, active.t0),
        active.death == PyVineyardMatrix::unpaired() ? inf : value_at(values0, values1, active.death, t1),
        active.birth,
        active.death,
        active.event0,
        event1
    });
}

void open_feature(VineyardLinearHomotopyResult& result,
                  ActiveVines& active_vines,
                  const Feature& feature,
                  double t0,
                  int event0)
{
    size_t vine_index = result.vines.size();
    result.vines.push_back(VineyardVine());
    active_vines[feature] = ActiveVine {
        vine_index,
        t0,
        event0,
        std::get<0>(feature),
        std::get<1>(feature)
    };
}

void reopen_feature(ActiveVines& active_vines,
                    const Feature& feature,
                    size_t vine_index,
                    double t0,
                    int event0)
{
    active_vines[feature] = ActiveVine {
        vine_index,
        t0,
        event0,
        std::get<0>(feature),
        std::get<1>(feature)
    };
}

ActiveVine close_feature(VineyardLinearHomotopyResult& result,
                         ActiveVines& active_vines,
                         const Feature& feature,
                         const std::vector<double>& values0,
                         const std::vector<double>& values1,
                         double t1,
                         int event1)
{
    auto it = active_vines.find(feature);
    if (it == active_vines.end())
        throw std::logic_error("vineyard linear homotopy active vine is missing");

    ActiveVine active = it->second;
    add_segment(result, active, values0, values1, t1, event1);
    active_vines.erase(it);
    return active;
}

void close_and_reopen_features(VineyardLinearHomotopyResult& result,
                               ActiveVines& active_vines,
                               const std::vector<Feature>& before,
                               const std::vector<Feature>& after,
                               const std::vector<double>& values0,
                               const std::vector<double>& values1,
                               double t,
                               int event)
{
    std::map<Feature, size_t> closed;
    for (const auto& feature : before)
    {
        ActiveVine active = close_feature(result, active_vines, feature, values0, values1, t, event);
        closed[feature] = active.vine_index;
    }
    for (const auto& feature : after)
    {
        auto it = closed.find(feature);
        if (it == closed.end())
            open_feature(result, active_vines, feature, t, event);
        else
            reopen_feature(active_vines, feature, it->second, t, event);
    }
}

void close_all_features(VineyardLinearHomotopyResult& result,
                        ActiveVines& active_vines,
                        const std::vector<double>& values0,
                        const std::vector<double>& values1,
                        double t,
                        int event)
{
    for (const auto& active : active_vines)
        add_segment(result, active.second, values0, values1, t, event);
    active_vines.clear();
}

bool boundary_contains(const Chain& column, Index row)
{
    auto it = std::lower_bound(column.begin(), column.end(), row,
                               [](const PyVineyardMatrix::Entry& entry, Index index)
                               { return entry.index() < index; });
    return it != column.end() && it->index() == row;
}

void validate_transposition(const Chains& boundary, Index first, Index second)
{
    if (boundary_contains(boundary[second], first))
        throw std::logic_error("vineyard linear homotopy transposition produced a non-filtration order");
}

template<class Vineyard>
void record_transposition(VineyardLinearHomotopyResult& result,
                          ActiveVines& active_vines,
                          Vineyard& vineyard,
                          const LinearHomotopyData& data,
                          double t,
                          Index position)
{
    Index first = vineyard.cell_at(position);
    Index second = vineyard.cell_at(position + 1);
    Index first_pair_before = vineyard.pair(first);
    Index second_pair_before = vineyard.pair(second);

    std::vector<Index> local_before;
    local_before.push_back(first);
    local_before.push_back(second);
    if (first_pair_before != Vineyard::unpaired())
        local_before.push_back(first_pair_before);
    if (second_pair_before != Vineyard::unpaired())
        local_before.push_back(second_pair_before);
    std::vector<Feature> before_features = local_features(vineyard, local_before);

    auto swapped = vineyard.transpose_position(position);
    Index first_pair_after = vineyard.pair(first);
    Index second_pair_after = vineyard.pair(second);

    std::vector<Index> local_after;
    local_after.push_back(first);
    local_after.push_back(second);
    if (first_pair_after != Vineyard::unpaired())
        local_after.push_back(first_pair_after);
    if (second_pair_after != Vineyard::unpaired())
        local_after.push_back(second_pair_after);
    std::vector<Feature> after_features = local_features(vineyard, local_after);

    int event = static_cast<int>(result.events.size());
    result.events.push_back(VineyardLinearHomotopyEvent {
        t,
        position,
        swapped.first,
        swapped.second,
        first_pair_before,
        second_pair_before,
        first_pair_after,
        second_pair_after
    });
    close_and_reopen_features(result, active_vines, before_features, after_features,
                              data.values0, data.values1, t, event);
}

template<class Vineyard>
void push_candidate(EventQueue& queue,
                    const Vineyard& vineyard,
                    const LinearHomotopyData& data,
                    double current_t,
                    Index position);

template<class Vineyard>
void push_candidate_neighborhood(EventQueue& queue,
                                 const Vineyard& vineyard,
                                 const LinearHomotopyData& data,
                                 double current_t,
                                 Index position)
{
    if (position > 0)
        push_candidate(queue, vineyard, data, current_t, position - 1);
    push_candidate(queue, vineyard, data, current_t, position);
    push_candidate(queue, vineyard, data, current_t, position + 1);
}

bool tie_less(const LinearHomotopyData& data, Index x, Index y)
{
    double sx = slope(data.values0, data.values1, x);
    double sy = slope(data.values0, data.values1, y);
    if (std::abs(sx - sy) > epsilon)
        return sx < sy;
    if (data.dimensions[x] != data.dimensions[y])
        return data.dimensions[x] < data.dimensions[y];
    return x < y;
}

bool adjacent_inverted_at(const LinearHomotopyData& data, Index left, Index right, double t)
{
    double left_value = value_at(data.values0, data.values1, left, t);
    double right_value = value_at(data.values0, data.values1, right, t);
    if (left_value > right_value + epsilon)
        return true;
    if (std::abs(left_value - right_value) <= epsilon)
        return tie_less(data, right, left);
    return false;
}

template<class Vineyard>
double crossing_time(const Vineyard& vineyard,
                     const LinearHomotopyData& data,
                     Index position,
                     double current_t)
{
    Index a = vineyard.cell_at(position);
    Index b = vineyard.cell_at(position + 1);
    double d = value_at(data.values0, data.values1, b, current_t) - value_at(data.values0, data.values1, a, current_t);
    double sd = slope(data.values0, data.values1, b) - slope(data.values0, data.values1, a);
    if (d > epsilon && sd < -epsilon)
    {
        double t = current_t - d / sd;
        if (t > current_t + epsilon && t <= 1.0 + epsilon)
            return std::min(1.0, t);
    }
    return std::numeric_limits<double>::infinity();
}

template<class Vineyard>
EventQueue build_event_queue(const Vineyard& vineyard,
                             const LinearHomotopyData& data,
                             double current_t)
{
    EventQueue queue;
    for (Index p = 0; p + 1 < vineyard.size(); ++p)
        push_candidate(queue, vineyard, data, current_t, p);
    return queue;
}

template<class Vineyard>
void push_candidate(EventQueue& queue,
                    const Vineyard& vineyard,
                    const LinearHomotopyData& data,
                    double current_t,
                    Index position)
{
    if (position + 1 >= vineyard.size())
        return;

    Index left = vineyard.cell_at(position);
    Index right = vineyard.cell_at(position + 1);
    double t = std::numeric_limits<double>::infinity();
    if (adjacent_inverted_at(data, left, right, current_t))
        t = current_t;
    else
        t = crossing_time(vineyard, data, position, current_t);

    if (std::isfinite(t))
        queue.push(CrossingCandidate {
            t,
            left,
            right,
            slope(data.values0, data.values1, right),
            data.dimensions[right],
            right
        });
}

template<class Vineyard>
bool pop_next_candidate(EventQueue& queue,
                        const Vineyard& vineyard,
                        const LinearHomotopyData& data,
                        double current_t,
                        CrossingCandidate& candidate,
                        Index& position)
{
    while (!queue.empty())
    {
        CrossingCandidate next = queue.top();
        queue.pop();

        if (next.time < current_t - epsilon)
            continue;
        if (next.time > 1.0 + epsilon)
            continue;
        if (vineyard.position(next.left) + 1 != vineyard.position(next.right))
            continue;
        double t = next.time <= current_t + epsilon ? current_t : std::min(1.0, next.time);
        position = vineyard.position(next.left);
        if (!adjacent_inverted_at(data, next.left, next.right, t))
            continue;
        candidate = next;
        candidate.time = t;
        return true;
    }
    return false;
}

template<class Vineyard>
VineyardLinearHomotopyResult run_linear_homotopy(const PyFiltration& filtration,
                                                const std::vector<double>& values0,
                                                const std::vector<double>& values1,
                                                const PyZpField& field)
{
    LinearHomotopyData data = prepare_linear_homotopy_data(filtration, values0, values1, field);
    auto* vineyard = new Vineyard(field, data.boundary);

    VineyardLinearHomotopyResult result;
    ActiveVines active_vines;
    for (const auto& feature : features(*vineyard))
        open_feature(result, active_vines, feature, 0.0, -1);

    double current_t = 0.0;
    EventQueue event_queue = build_event_queue(*vineyard, data, current_t);

    while (current_t < 1.0 - epsilon)
    {
        CrossingCandidate candidate;
        Index position = 0;
        if (!pop_next_candidate(event_queue, *vineyard, data, current_t, candidate, position) || candidate.time > 1.0 - epsilon)
        {
            close_all_features(result, active_vines, data.values0, data.values1, 1.0, -1);
            current_t = 1.0;
            break;
        }

        current_t = candidate.time;
        Index first = vineyard->cell_at(position);
        Index second = vineyard->cell_at(position + 1);
        validate_transposition(data.boundary, first, second);
        record_transposition(result, active_vines, *vineyard, data, current_t, position);
        push_candidate_neighborhood(event_queue, *vineyard, data, current_t, position);
    }

    result.final_order = current_order(*vineyard);
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

    py::class_<PyVineyardV>(m, "VineyardV", "matrix_v-based vineyard state for adjacent transpositions")
        .def(py::init([make_chains](Columns columns, PyZpField field)
                        {
                            return new PyVineyardV(field, make_chains(std::move(columns)));
                        }), "columns"_a = Columns(), "field"_a = PyZpField(2))
        .def(py::init([](const PyFiltration& filtration, PyZpField field)
                       {
                           return new PyVineyardV(field, boundary_from_filtration(filtration, field));
                       }), "filtration"_a, "field"_a = PyZpField(2))
        .def("__len__",      &PyVineyardV::size,      "number of cells")
        .def("field",        &PyVineyardV::field,     "field used for coefficients")
        .def("cell_at",      &PyVineyardV::cell_at,   "cell id at a current filtration position")
        .def("position",     &PyVineyardV::position,  "current filtration position of a cell id")
        .def("low",          &PyVineyardV::low,       "cached low cell id of a reduced column")
        .def("pivot",        &PyVineyardV::pivot,     "column id with the given low cell id")
        .def("pair",         &PyVineyardV::pair,      "persistence pair of the given cell id")
        .def("transpose_position", &PyVineyardV::transpose_position, "position"_a,
                                                             "repair the vineyard state after transposing adjacent filtration positions")
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
                           return new PyVineyardU(field, boundary_from_filtration(filtration, field));
                       }), "filtration"_a, "field"_a = PyZpField(2))
        .def("__len__",      &PyVineyardU::size,      "number of cells")
        .def("field",        &PyVineyardU::field,     "field used for coefficients")
        .def("cell_at",      &PyVineyardU::cell_at,   "cell id at a current filtration position")
        .def("position",     &PyVineyardU::position,  "current filtration position of a cell id")
        .def("low",          &PyVineyardU::low,       "cached low cell id of a reduced column")
        .def("pivot",        &PyVineyardU::pivot,     "column id with the given low cell id")
        .def("pair",         &PyVineyardU::pair,      "persistence pair of the given cell id")
        .def("transpose_position", &PyVineyardU::transpose_position, "position"_a,
                                                             "repair the vineyard state after transposing adjacent filtration positions")
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
                  return py::cast(new PyVineyardV(field, boundary_from_filtration(filtration, field)), py::return_value_policy::take_ownership);
              if (method == "matrix_u")
                  return py::cast(new PyVineyardU(field, boundary_from_filtration(filtration, field)), py::return_value_policy::take_ownership);
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
          "compute a vineyard linear homotopy between two filtration functions");
}
