#ifndef DIONYSUS_FREUDENTHAL_H
#define DIONYSUS_FREUDENTHAL_H

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <vector>

#include "rips.h"
#include "simplex.h"

namespace dionysus
{

namespace detail
{

struct FreudenthalDummyDistances
{
    using IndexType = int;
    using DistanceType = float;

    DistanceType operator()(int, int) const { return 0; }

    IndexType begin() const { return 0; }
    IndexType end() const { return 0; }
};

inline std::vector<Simplex<std::vector<int>>> freudenthal_delta_simplices(size_t dimension)
{
    using Delta = std::vector<int>;
    using DeltaSimplex = Simplex<Delta>;
    using Rips = dionysus::Rips<FreudenthalDummyDistances, DeltaSimplex>;

    Rips::VertexContainer vertices;
    vertices.emplace_back(dimension, 0);

    Rips::VertexContainer candidates;
    for (size_t i = 1; i < (1u << dimension); ++i)
    {
        candidates.emplace_back();
        for (size_t j = 0; j < dimension; ++j)
            candidates.back().emplace_back((i >> j) & 1 ? 1 : 0);
    }

    std::vector<DeltaSimplex> delta_simplices;
    Rips::bron_kerbosch(vertices, candidates, std::prev(candidates.begin()), dimension,
                        [](const Delta& v1, const Delta& v2)
                        {
                            bool pos = false, neg = false;
                            for (size_t i = 0; i < v1.size(); ++i)
                            {
                                auto diff = v1[i] - v2[i];
                                if (diff < 0)
                                    neg = true;
                                else if (diff > 0)
                                    pos = true;
                            }

                            return (pos ^ neg);
                        },
                        [&delta_simplices](const DeltaSimplex& s)
                        {
                            delta_simplices.push_back(s);
                        });
    return delta_simplices;
}

}

template<class Filtration, class Value, class Shape, class Strides, class Compare>
Filtration fill_freudenthal(const Shape& shape, const Strides& strides, const Value* data, bool reverse, Compare cmp)
{
    using Delta = std::vector<int>;
    using Cell = typename Filtration::Cell;
    using CellData = typename Cell::Data;

    Filtration filtration;
    auto delta_simplices = detail::freudenthal_delta_simplices(shape.size());

    Delta v(shape.size(), 0);
    while (v.back() <= shape.back())
    {
        for (const auto& s : delta_simplices)
        {
            std::vector<size_t> vertices;
            for (auto& u : s)
            {
                size_t uidx = 0;
                size_t i = 0;
                for (; i < v.size(); ++i)
                {
                    auto nbr = v[i] + u[i];
                    if (nbr >= shape[i])
                        break;
                    uidx += nbr * strides[i];
                }
                if (i == v.size())
                    vertices.push_back(uidx);
            }

            if (vertices.size() == s.dimension() + 1)
            {
                auto selected = *std::max_element(vertices.begin(), vertices.end(),
                                                  [data, reverse](size_t u, size_t v)
                                                  {
                                                      return reverse ? data[v] < data[u] : data[u] < data[v];
                                                  });
                filtration.emplace_back(vertices, static_cast<CellData>(data[selected]));
            }
        }

        size_t i = 0;
        while(i < v.size())
        {
            if (++v[i] < shape[i] || i == v.size() - 1)
                break;

            v[i] = 0;
            ++i;
        }
    }

    filtration.sort(cmp);
    return filtration;
}

}

#endif
