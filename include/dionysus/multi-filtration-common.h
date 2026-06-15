#ifndef DIONYSUS_MULTI_FILTRATION_COMMON_H
#define DIONYSUS_MULTI_FILTRATION_COMMON_H

#include <cstddef>

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

}

}

#endif
