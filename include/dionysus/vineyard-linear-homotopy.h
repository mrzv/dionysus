#ifndef DIONYSUS_VINEYARD_LINEAR_HOMOTOPY_H
#define DIONYSUS_VINEYARD_LINEAR_HOMOTOPY_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace dionysus
{

constexpr double vineyard_linear_homotopy_epsilon = 1e-10;

template<class Index>
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
    bool    pairing_switched;
};

template<class Index>
struct VineyardLinearSegment
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

template<class Index>
struct VineyardLinearVine
{
    std::vector<VineyardLinearSegment<Index>> segments;
};

template<class Index>
struct VineyardLinearActiveVine
{
    size_t  vine_index;
    double  t0;
    int     event0;
    Index   birth;
    Index   death;
};

template<class Chains>
struct VineyardLinearHomotopyData
{
    using Index = typename Chains::value_type::value_type::Index;

    Chains                  boundary;
    std::vector<Index>      stable_to_input;
    std::vector<Index>      input_to_stable;
    std::vector<double>     values0;
    std::vector<double>     values1;
    std::vector<unsigned>   dimensions;
};

template<class Index>
struct VineyardLinearCrossingCandidate
{
    double      time;
    Index       left;
    Index       right;
    double      priority_slope;
    unsigned    priority_dimension;
    Index       priority_cell;

    bool operator<(const VineyardLinearCrossingCandidate& other) const
    {
        if (time != other.time)
            return time > other.time;
        if (priority_slope != other.priority_slope)
            return priority_slope > other.priority_slope;
        if (priority_dimension != other.priority_dimension)
            return priority_dimension > other.priority_dimension;
        return priority_cell > other.priority_cell;
    }
};

template<class Vineyard>
std::vector<typename Vineyard::Index> current_vineyard_linear_order(const Vineyard& vineyard)
{
    using Index = typename Vineyard::Index;

    std::vector<Index> order;
    order.reserve(vineyard.size());
    for (Index p = 0; p < vineyard.size(); ++p)
        order.push_back(vineyard.cell_at(p));
    return order;
}

template<class Vineyard, class Feature>
bool vineyard_linear_feature_for_cell(const Vineyard& vineyard, typename Vineyard::Index cell, Feature& feature)
{
    if (cell == Vineyard::unpaired())
        return false;

    auto pair = vineyard.pair(cell);
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

template<class Vineyard, class Feature>
std::vector<Feature> vineyard_linear_features(const Vineyard& vineyard)
{
    using Index = typename Vineyard::Index;

    std::vector<Feature> result;
    for (Index cell = 0; cell < vineyard.size(); ++cell)
    {
        Feature feature;
        if (vineyard_linear_feature_for_cell(vineyard, cell, feature))
            result.push_back(feature);
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

template<class Vineyard, class Feature>
std::vector<Feature> local_vineyard_linear_features(const Vineyard& vineyard, const std::vector<typename Vineyard::Index>& cells)
{
    std::vector<Feature> result;
    for (auto cell : cells)
    {
        Feature feature;
        if (vineyard_linear_feature_for_cell(vineyard, cell, feature))
            result.push_back(feature);
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

inline double vineyard_linear_value_at(const std::vector<double>& values0,
                                       const std::vector<double>& values1,
                                       size_t cell,
                                       double t)
{
    return values0[cell] + t * (values1[cell] - values0[cell]);
}

inline double vineyard_linear_slope(const std::vector<double>& values0,
                                    const std::vector<double>& values1,
                                    size_t cell)
{
    return values1[cell] - values0[cell];
}

template<class Result, class ActiveVine, class Index>
void add_vineyard_linear_segment(Result& result,
                                 const ActiveVine& active,
                                 const std::vector<double>& values0,
                                 const std::vector<double>& values1,
                                 double t1,
                                 int event1,
                                 Index unpaired)
{
    if (t1 <= active.t0)
        return;

    double inf = std::numeric_limits<double>::infinity();
    result.vines[active.vine_index].segments.push_back({
        active.t0,
        t1,
        vineyard_linear_value_at(values0, values1, active.birth, active.t0),
        vineyard_linear_value_at(values0, values1, active.birth, t1),
        active.death == unpaired ? inf : vineyard_linear_value_at(values0, values1, active.death, active.t0),
        active.death == unpaired ? inf : vineyard_linear_value_at(values0, values1, active.death, t1),
        active.birth,
        active.death,
        active.event0,
        event1
    });
}

template<class Result, class ActiveVines, class Feature>
void open_vineyard_linear_feature(Result& result,
                                  ActiveVines& active_vines,
                                  const Feature& feature,
                                  double t0,
                                  int event0)
{
    size_t vine_index = result.vines.size();
    result.vines.emplace_back();
    active_vines[feature] = typename ActiveVines::mapped_type {
        vine_index,
        t0,
        event0,
        std::get<0>(feature),
        std::get<1>(feature)
    };
}

template<class ActiveVines, class Feature>
void reopen_vineyard_linear_feature(ActiveVines& active_vines,
                                    const Feature& feature,
                                    size_t vine_index,
                                    double t0,
                                    int event0)
{
    active_vines[feature] = typename ActiveVines::mapped_type {
        vine_index,
        t0,
        event0,
        std::get<0>(feature),
        std::get<1>(feature)
    };
}

template<class Result, class ActiveVines, class Feature, class Index>
typename ActiveVines::mapped_type close_vineyard_linear_feature(Result& result,
                                                                ActiveVines& active_vines,
                                                                const Feature& feature,
                                                                const std::vector<double>& values0,
                                                                const std::vector<double>& values1,
                                                                double t1,
                                                                int event1,
                                                                Index unpaired)
{
    auto it = active_vines.find(feature);
    if (it == active_vines.end())
        throw std::logic_error("vineyard linear homotopy active vine is missing");

    auto active = it->second;
    add_vineyard_linear_segment(result, active, values0, values1, t1, event1, unpaired);
    active_vines.erase(it);
    return active;
}

template<class Feature, class Index>
Feature transpose_vineyard_linear_feature(const Feature& feature, Index first, Index second)
{
    auto transpose_cell = [first, second](Index cell)
    {
        if (cell == first)
            return second;
        if (cell == second)
            return first;
        return cell;
    };

    return std::make_tuple(transpose_cell(std::get<0>(feature)), transpose_cell(std::get<1>(feature)));
}

template<class Feature>
bool contains_vineyard_linear_feature(const std::vector<Feature>& features, const Feature& feature)
{
    return std::binary_search(features.begin(), features.end(), feature);
}

template<class Result, class ActiveVines, class Feature, class Index>
void continue_vineyard_linear_features(Result& result,
                                       ActiveVines& active_vines,
                                       const std::vector<Feature>& before,
                                       const std::vector<Feature>& after,
                                       const std::vector<double>& values0,
                                       const std::vector<double>& values1,
                                       double t,
                                       int event,
                                       Index first,
                                       Index second,
                                       bool pairing_switched,
                                       Index unpaired)
{
    using ActiveVine = typename ActiveVines::mapped_type;
    using Continuation = std::pair<Feature, ActiveVine>;

    std::vector<Feature> continued;
    std::vector<Continuation> reopen;
    continued.reserve(before.size());
    reopen.reserve(before.size());

    for (const auto& feature : before)
    {
        Feature target = pairing_switched ? transpose_vineyard_linear_feature(feature, first, second) : feature;
        if (!contains_vineyard_linear_feature(after, target))
        {
            close_vineyard_linear_feature(result, active_vines, feature, values0, values1, t, event, unpaired);
            continue;
        }

        continued.push_back(target);
        if (target == feature)
            continue;

        auto active = close_vineyard_linear_feature(result, active_vines, feature, values0, values1, t, event, unpaired);
        reopen.emplace_back(target, active);
    }

    for (const auto& continuation : reopen)
    {
        if (active_vines.find(continuation.first) != active_vines.end())
            throw std::logic_error("vineyard linear homotopy feature continuation collision");
        reopen_vineyard_linear_feature(active_vines, continuation.first, continuation.second.vine_index, t, event);
    }

    std::sort(continued.begin(), continued.end());
    continued.erase(std::unique(continued.begin(), continued.end()), continued.end());
    for (const auto& feature : after)
        if (!contains_vineyard_linear_feature(continued, feature))
            open_vineyard_linear_feature(result, active_vines, feature, t, event);
}

template<class Result, class ActiveVines, class Index>
void close_all_vineyard_linear_features(Result& result,
                                        ActiveVines& active_vines,
                                        const std::vector<double>& values0,
                                        const std::vector<double>& values1,
                                        double t,
                                        int event,
                                        Index unpaired)
{
    for (const auto& active : active_vines)
        add_vineyard_linear_segment(result, active.second, values0, values1, t, event, unpaired);
    active_vines.clear();
}

template<class Chain, class Index>
bool vineyard_linear_boundary_contains(const Chain& column, Index row)
{
    auto it = std::lower_bound(column.begin(), column.end(), row,
                               [](const typename Chain::value_type& entry, Index index)
                               { return entry.index() < index; });
    return it != column.end() && it->index() == row;
}

template<class Chains, class Index>
bool vineyard_linear_transposition_preserves_filtration(const Chains& boundary, Index first, Index second)
{
    return !vineyard_linear_boundary_contains(boundary[second], first);
}

template<class Chains, class Index>
void validate_vineyard_linear_transposition(const Chains& boundary, Index first, Index second)
{
    if (!vineyard_linear_transposition_preserves_filtration(boundary, first, second))
        throw std::logic_error("vineyard linear homotopy transposition produced a non-filtration order");
}

template<class Result, class ActiveVines, class Vineyard, class Data, class Feature>
void record_vineyard_linear_transposition(Result& result,
                                          ActiveVines& active_vines,
                                          Vineyard& vineyard,
                                          const Data& data,
                                          double t,
                                          typename Data::Index position)
{
    using Index = typename Data::Index;

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
    auto before_features = local_vineyard_linear_features<Vineyard, Feature>(vineyard, local_before);

    auto swapped = vineyard.transpose(position);
    bool pairing_switched = std::get<2>(swapped);
    Index first_pair_after = vineyard.pair(first);
    Index second_pair_after = vineyard.pair(second);

    std::vector<Index> local_after;
    local_after.push_back(first);
    local_after.push_back(second);
    if (first_pair_after != Vineyard::unpaired())
        local_after.push_back(first_pair_after);
    if (second_pair_after != Vineyard::unpaired())
        local_after.push_back(second_pair_after);
    auto after_features = local_vineyard_linear_features<Vineyard, Feature>(vineyard, local_after);

    int event = static_cast<int>(result.events.size());
    result.events.push_back({
        t,
        position,
        std::get<0>(swapped),
        std::get<1>(swapped),
        first_pair_before,
        second_pair_before,
        first_pair_after,
        second_pair_after,
        pairing_switched
    });
    continue_vineyard_linear_features<Result, ActiveVines, Feature>(result, active_vines, before_features, after_features,
                                                                    data.values0, data.values1, t, event, first, second,
                                                                    pairing_switched, Vineyard::unpaired());
}

template<class Data>
bool vineyard_linear_tie_less(const Data& data, typename Data::Index x, typename Data::Index y)
{
    double sx = vineyard_linear_slope(data.values0, data.values1, x);
    double sy = vineyard_linear_slope(data.values0, data.values1, y);
    if (sx != sy)
        return sx < sy;
    if (data.dimensions[x] != data.dimensions[y])
        return data.dimensions[x] < data.dimensions[y];
    return x < y;
}

template<class Data>
double vineyard_linear_crossing_time(const Data& data,
                                     typename Data::Index left,
                                     typename Data::Index right,
                                     double current_t)
{
    double d = vineyard_linear_value_at(data.values0, data.values1, right, current_t) -
               vineyard_linear_value_at(data.values0, data.values1, left, current_t);
    double sd = vineyard_linear_slope(data.values0, data.values1, right) -
                vineyard_linear_slope(data.values0, data.values1, left);
    if (d >= 0.0 && sd < 0.0)
    {
        double t = current_t - d / sd;
        if (t > current_t && t <= 1.0)
            return t;
    }
    return std::numeric_limits<double>::infinity();
}

template<class Vineyard, class Data>
bool vineyard_linear_candidate_time(const Vineyard& vineyard,
                                    const Data& data,
                                    typename Data::Index position,
                                    double current_t,
                                    double& candidate_t)
{
    auto left = vineyard.cell_at(position);
    auto right = vineyard.cell_at(position + 1);
    double left_value = vineyard_linear_value_at(data.values0, data.values1, left, current_t);
    double right_value = vineyard_linear_value_at(data.values0, data.values1, right, current_t);
    if (left_value > right_value || (left_value == right_value && vineyard_linear_tie_less(data, right, left)))
    {
        candidate_t = current_t;
        return true;
    }

    double crossing = vineyard_linear_crossing_time(data, left, right, current_t);
    if (!std::isfinite(crossing))
        return false;

    candidate_t = crossing;
    return true;
}

template<class EventQueue, class Vineyard, class Data>
void push_vineyard_linear_candidate(EventQueue& queue,
                                     const Vineyard& vineyard,
                                     const Data& data,
                                     double current_t,
                                     typename Data::Index position)
{
    if (position + 1 >= vineyard.size())
        return;

    auto left = vineyard.cell_at(position);
    auto right = vineyard.cell_at(position + 1);
    double time = std::numeric_limits<double>::infinity();
    if (vineyard_linear_candidate_time(vineyard, data, position, current_t, time))
        queue.push(typename EventQueue::value_type {
            time,
            left,
            right,
            vineyard_linear_slope(data.values0, data.values1, right),
            data.dimensions[right],
            right
        });
}

template<class EventQueue, class Vineyard, class Data>
void push_vineyard_linear_candidate_neighborhood(EventQueue& queue,
                                                  const Vineyard& vineyard,
                                                  const Data& data,
                                                  double current_t,
                                                  typename Data::Index position)
{
    if (position > 0)
        push_vineyard_linear_candidate(queue, vineyard, data, current_t, position - 1);
    push_vineyard_linear_candidate(queue, vineyard, data, current_t, position);
    push_vineyard_linear_candidate(queue, vineyard, data, current_t, position + 1);
}

template<class EventQueue, class Vineyard, class Data>
EventQueue build_vineyard_linear_event_queue(const Vineyard& vineyard,
                                              const Data& data,
                                              double current_t)
{
    EventQueue queue;
    for (typename Data::Index p = 0; p + 1 < vineyard.size(); ++p)
        push_vineyard_linear_candidate(queue, vineyard, data, current_t, p);
    return queue;
}

template<class EventQueue, class Vineyard, class Data>
bool pop_next_vineyard_linear_candidate(EventQueue& queue,
                                         const Vineyard& vineyard,
                                         const Data& data,
                                         double current_t,
                                         typename EventQueue::value_type& candidate,
                                         typename Data::Index& position)
{
    while (!queue.empty())
    {
        auto next = queue.top();
        queue.pop();

        if (next.time < current_t)
            continue;
        if (next.time > 1.0)
            continue;
        if (vineyard.position(next.left) + 1 != vineyard.position(next.right))
            continue;
        if (!vineyard_linear_transposition_preserves_filtration(data.boundary, next.left, next.right))
            continue;
        position = vineyard.position(next.left);
        double recomputed_time = std::numeric_limits<double>::infinity();
        if (!vineyard_linear_candidate_time(vineyard, data, position, current_t, recomputed_time))
            continue;
        if (recomputed_time != next.time)
        {
            queue.push(typename EventQueue::value_type {
                recomputed_time,
                next.left,
                next.right,
                vineyard_linear_slope(data.values0, data.values1, next.right),
                data.dimensions[next.right],
                next.right
            });
            continue;
        }
        candidate = next;
        candidate.time = recomputed_time;
        return true;
    }
    return false;
}

template<class Filtration>
bool vineyard_linear_filtration_less(const Filtration& filtration,
                                     const std::vector<double>& values,
                                     size_t x,
                                     size_t y)
{
    if (values[x] != values[y])
        return values[x] < values[y];
    return filtration[x] < filtration[y];
}

template<class Filtration>
void validate_vineyard_linear_endpoint(const Filtration& filtration,
                                       const std::vector<double>& values,
                                       double epsilon = vineyard_linear_homotopy_epsilon)
{
    for (size_t i = 0; i < filtration.size(); ++i)
        for (auto it = filtration[i].boundary_begin(); it != filtration[i].boundary_end(); ++it)
        {
            size_t face = filtration.index(*it, filtration.size());
            if (values[face] > values[i] + epsilon)
                throw std::invalid_argument("linear homotopy values are not a filtration on the complex");
        }
}

template<class Chains, class Filtration, class Field>
VineyardLinearHomotopyData<Chains> prepare_vineyard_linear_homotopy_data(const Filtration& filtration,
                                                                         const std::vector<double>& input_values0,
                                                                         const std::vector<double>& input_values1,
                                                                         const Field& field)
{
    using Data = VineyardLinearHomotopyData<Chains>;
    using Chain = typename Chains::value_type;
    using Entry = typename Chain::value_type;
    using Index = typename Entry::Index;

    if (input_values0.size() != filtration.size() || input_values1.size() != filtration.size())
        throw std::invalid_argument("linear homotopy value vectors must match filtration size");

    validate_vineyard_linear_endpoint(filtration, input_values0);
    validate_vineyard_linear_endpoint(filtration, input_values1);

    Data data;
    data.stable_to_input.resize(filtration.size());
    data.input_to_stable.resize(filtration.size());
    for (Index i = 0; i < filtration.size(); ++i)
        data.stable_to_input[i] = i;

    std::sort(data.stable_to_input.begin(), data.stable_to_input.end(),
              [&filtration, &input_values0](Index x, Index y)
              { return vineyard_linear_filtration_less(filtration, input_values0, x, y); });

    for (Index stable = 0; stable < data.stable_to_input.size(); ++stable)
        data.input_to_stable[data.stable_to_input[stable]] = stable;

    data.values0.resize(filtration.size());
    data.values1.resize(filtration.size());
    data.dimensions.resize(filtration.size());
    data.boundary.resize(filtration.size());

    for (Index stable = 0; stable < data.stable_to_input.size(); ++stable)
    {
        Index input = data.stable_to_input[stable];
        const auto& simplex = filtration[input];
        data.values0[stable] = input_values0[input];
        data.values1[stable] = input_values1[input];
        data.dimensions[stable] = simplex.dimension();

        Chain column;
        for (auto it = simplex.boundary_begin(field); it != simplex.boundary_end(field); ++it)
        {
            Index face_input = filtration.index(it->index(), filtration.size());
            column.emplace_back(it->element(), data.input_to_stable[face_input]);
        }

        std::sort(column.begin(), column.end(), [](const Entry& x, const Entry& y)
                  { return x.index() < y.index(); });
        data.boundary[stable] = std::move(column);
    }

    return data;
}

}

#endif
