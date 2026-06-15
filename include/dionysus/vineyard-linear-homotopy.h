#ifndef DIONYSUS_VINEYARD_LINEAR_HOMOTOPY_H
#define DIONYSUS_VINEYARD_LINEAR_HOMOTOPY_H

#include <algorithm>
#include <cmath>
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
