#ifndef DIONYSUS_MULTI_FILTRATION_COMMON_H
#define DIONYSUS_MULTI_FILTRATION_COMMON_H

#include <cstddef>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace dionysus
{

namespace detail
{

template<class Cell>
bool indexed_cell_less(const Cell& c, size_t i, const Cell& other, size_t other_i)
{
    return c < other || (c == other && i < other_i);
}

template<class IndexedCell>
struct CellWithoutIndex
{
    using Cell = typename IndexedCell::Cell;

    Cell&       operator()(IndexedCell& c) const       { return c; }
    const Cell& operator()(const IndexedCell& c) const { return c; }
};

template<class Order>
std::vector<std::reference_wrapper<const typename Order::value_type>>
rearrangement_references(const Order& order, const std::vector<size_t>& indices)
{
    std::vector<std::reference_wrapper<const typename Order::value_type>> references;
    references.reserve(indices.size());
    for (size_t i : indices)
        references.push_back(std::cref(order[i]));
    return references;
}

inline std::vector<size_t> reverse_rearrangement_indices(const std::vector<size_t>& indices)
{
    std::vector<size_t> reverse_indices(indices.size());
    size_t i = 0;
    for (size_t idx : indices)
        reverse_indices[idx] = i++;
    return reverse_indices;
}

template<class Order>
std::vector<size_t> current_to_sorted_indices(const Order& order, size_t size)
{
    size_t i = 0;
    std::vector<size_t> indices(size);
    for (auto& c : order)
        indices[c.i] = i++;
    return indices;
}

template<class Order, class Complex, class ProjectComplex, class UpdateCell>
void update_order_indices(Order& order, Complex& complex, ProjectComplex project_complex, UpdateCell update_cell)
{
    size_t i = 0;
    for(auto it = order.begin(); it != order.end(); ++it)
    {
        auto cit = project_complex(it);
        complex.modify(cit, [i,&update_cell](typename Order::value_type& c) { update_cell(c, i); });
        ++i;
    }
}

template<class Complex, class IndexedCell>
typename Complex::const_iterator preceding_indexed_cell(const Complex& complex, const IndexedCell& query)
{
    auto it = complex.upper_bound(query);
    if (it == complex.begin())
        return complex.end();
    --it;
    return it;
}

template<bool checked_index, class Iterator, class Cell>
void validate_indexed_cell_lookup(Iterator candidate, Iterator end, const Cell& cell, size_t index)
{
    if (checked_index && (candidate == end || !(*candidate == cell)))
    {
        std::ostringstream oss;
        oss << "Trying to access non-existent cell: " << cell << ' ' << index;
        throw std::runtime_error(oss.str());
    }
}

}

}

#endif
