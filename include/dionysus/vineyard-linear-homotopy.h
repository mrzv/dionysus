#ifndef DIONYSUS_VINEYARD_LINEAR_HOMOTOPY_H
#define DIONYSUS_VINEYARD_LINEAR_HOMOTOPY_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace dionysus
{

constexpr double vineyard_linear_homotopy_epsilon = 1e-10;

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
        if (std::abs(time - other.time) > vineyard_linear_homotopy_epsilon)
            return time > other.time;
        if (std::abs(priority_slope - other.priority_slope) > vineyard_linear_homotopy_epsilon)
            return priority_slope > other.priority_slope;
        if (priority_dimension != other.priority_dimension)
            return priority_dimension > other.priority_dimension;
        return priority_cell > other.priority_cell;
    }
};

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

template<class Data>
bool vineyard_linear_tie_less(const Data& data, typename Data::Index x, typename Data::Index y)
{
    double sx = vineyard_linear_slope(data.values0, data.values1, x);
    double sy = vineyard_linear_slope(data.values0, data.values1, y);
    if (std::abs(sx - sy) > vineyard_linear_homotopy_epsilon)
        return sx < sy;
    if (data.dimensions[x] != data.dimensions[y])
        return data.dimensions[x] < data.dimensions[y];
    return x < y;
}

template<class Data>
bool vineyard_linear_adjacent_inverted_at(const Data& data,
                                         typename Data::Index left,
                                         typename Data::Index right,
                                         double t)
{
    double left_value = vineyard_linear_value_at(data.values0, data.values1, left, t);
    double right_value = vineyard_linear_value_at(data.values0, data.values1, right, t);
    if (left_value > right_value + vineyard_linear_homotopy_epsilon)
        return true;
    if (std::abs(left_value - right_value) <= vineyard_linear_homotopy_epsilon)
        return vineyard_linear_tie_less(data, right, left);
    return false;
}

template<class Vineyard, class Data>
double vineyard_linear_crossing_time(const Vineyard& vineyard,
                                     const Data& data,
                                     typename Data::Index position,
                                     double current_t)
{
    auto a = vineyard.cell_at(position);
    auto b = vineyard.cell_at(position + 1);
    double d = vineyard_linear_value_at(data.values0, data.values1, b, current_t) - vineyard_linear_value_at(data.values0, data.values1, a, current_t);
    double sd = vineyard_linear_slope(data.values0, data.values1, b) - vineyard_linear_slope(data.values0, data.values1, a);
    if (d > vineyard_linear_homotopy_epsilon && sd < -vineyard_linear_homotopy_epsilon)
    {
        double t = current_t - d / sd;
        if (t > current_t + vineyard_linear_homotopy_epsilon && t <= 1.0 + vineyard_linear_homotopy_epsilon)
            return std::min(1.0, t);
    }
    return std::numeric_limits<double>::infinity();
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
    double t = std::numeric_limits<double>::infinity();
    if (vineyard_linear_adjacent_inverted_at(data, left, right, current_t))
        t = current_t;
    else
        t = vineyard_linear_crossing_time(vineyard, data, position, current_t);

    if (std::isfinite(t))
        queue.push(typename EventQueue::value_type {
            t,
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

        if (next.time < current_t - vineyard_linear_homotopy_epsilon)
            continue;
        if (next.time > 1.0 + vineyard_linear_homotopy_epsilon)
            continue;
        if (vineyard.position(next.left) + 1 != vineyard.position(next.right))
            continue;
        double t = next.time <= current_t + vineyard_linear_homotopy_epsilon ? current_t : std::min(1.0, next.time);
        position = vineyard.position(next.left);
        if (!vineyard_linear_adjacent_inverted_at(data, next.left, next.right, t))
            continue;
        candidate = next;
        candidate.time = t;
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
