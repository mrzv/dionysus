#ifndef DIONYSUS_MULTI_FILTRATION_COMMON_H
#define DIONYSUS_MULTI_FILTRATION_COMMON_H

#include <cstddef>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/random_access_index.hpp>

namespace dionysus
{

namespace detail
{

template<class Cell>
bool indexed_cell_less(const Cell& c, size_t i, const Cell& other, size_t other_i);

template<class Cell_>
struct DuplicateCellWithIndex: Cell_
{
    using Cell = Cell_;

    DuplicateCellWithIndex(const Cell& c_, size_t i_):
        Cell(c_), i(i_)             {}
    DuplicateCellWithIndex(Cell&& c_, size_t i_):
        Cell(std::move(c_)), i(i_)  {}

    size_t  i;

    friend bool operator<(const DuplicateCellWithIndex& c, const DuplicateCellWithIndex& other)
    {
        return indexed_cell_less(static_cast<const Cell&>(c), c.i,
                                 static_cast<const Cell&>(other), other.i);
    }
};

template<class Cell_>
struct LinkedDuplicateCellWithIndex: Cell_
{
    using Cell = Cell_;

    LinkedDuplicateCellWithIndex(const Cell& c_, size_t i_, size_t l_ = 0):
        Cell(c_), i(i_), linked(l_)             {}
    LinkedDuplicateCellWithIndex(Cell&& c_, size_t i_, size_t l_ = 0):
        Cell(std::move(c_)), i(i_), linked(l_)  {}

    size_t  i;
    size_t  linked;

    friend bool operator<(const LinkedDuplicateCellWithIndex& c, const LinkedDuplicateCellWithIndex& other)
    {
        return indexed_cell_less(static_cast<const Cell&>(c), c.i,
                                 static_cast<const Cell&>(other), other.i);
    }
};

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

template<class IndexedCell, bool checked_index>
class DuplicateFiltrationStorage
{
    public:
        struct order {};

        using Cell = typename IndexedCell::Cell;
        using CellLookupIndex = boost::multi_index::ordered_non_unique<boost::multi_index::identity<IndexedCell>>;
        using Container = boost::multi_index_container<IndexedCell,
                                                       boost::multi_index::indexed_by<CellLookupIndex,
                                                                                       boost::multi_index::random_access<boost::multi_index::tag<order>>>>;
        using value_type = typename Container::value_type;

        using Complex = typename Container::template nth_index<0>::type;
        using Order = typename Container::template nth_index<1>::type;

        using OrderCell = CellWithoutIndex<IndexedCell>;
        using OrderConstIterator = boost::transform_iterator<OrderCell, typename Order::const_iterator>;
        using OrderIterator = boost::transform_iterator<OrderCell, typename Order::iterator>;

        DuplicateFiltrationStorage() = default;
        DuplicateFiltrationStorage(DuplicateFiltrationStorage&& other) = default;
        DuplicateFiltrationStorage& operator=(DuplicateFiltrationStorage&& other) = default;

        const Cell&         operator[](size_t i) const          { return get_order()[i]; }
        bool                contains(const Cell& s) const       { return get_complex().contains(s, std::less<Cell>()); }

        OrderConstIterator  begin() const                       { return OrderConstIterator(get_order().begin()); }
        OrderConstIterator  end() const                         { return OrderConstIterator(get_order().end()); }
        OrderIterator       begin()                             { return OrderIterator(get_order().begin()); }
        OrderIterator       end()                               { return OrderIterator(get_order().end()); }
        size_t              size() const                        { return cells_.size(); }
        void                clear()                             { return Container().swap(cells_); }

        Cell&               back()                              { return const_cast<Cell&>(get_order().back()); }
        const Cell&         back() const                        { return get_order().back(); }

    protected:
        const Complex&      get_complex() const                 { return cells_.template get<0>(); }
        Complex&            get_complex()                       { return cells_.template get<0>(); }
        const Order&        get_order() const                   { return cells_.template get<order>(); }
        Order&              get_order()                         { return cells_.template get<order>(); }
        size_t              next_index() const                  { return cells_.size(); }

        template<class Iterator>
        typename Order::const_iterator
                            project_order(Iterator it) const    { return cells_.template project<order>(it); }
        template<class Iterator>
        typename Complex::iterator
                            project_complex(Iterator it)        { return cells_.template project<0>(it); }

        size_t preceding_index(const Cell& s, size_t i) const
        {
            auto it = preceding_indexed_cell(get_complex(), IndexedCell(s, i));
            validate_indexed_cell_lookup<checked_index>(it, get_complex().end(), s, i);
            return project_order(it) - get_order().begin();
        }

        void push_indexed(const IndexedCell& c)                  { get_order().push_back(c); }
        void push_indexed(IndexedCell&& c)                       { get_order().push_back(std::move(c)); }
        void replace_indexed(size_t i, const IndexedCell& c)     { get_order().replace(get_order().begin() + i, c); }

        template<class UpdateCell>
        void rearrange_with_updates(const std::vector<size_t>& indices, UpdateCell update_cell)
        {
            auto references = rearrangement_references(get_order(), indices);
            get_order().rearrange(references.begin());
            update_indices(update_cell);
        }

        template<class UpdateCell>
        void rearrange_with_index_map(const std::vector<size_t>& indices, UpdateCell update_cell)
        {
            auto references = rearrangement_references(get_order(), indices);
            auto reverse_indices = reverse_rearrangement_indices(indices);
            get_order().rearrange(references.begin());
            update_indices(reverse_indices, update_cell);
        }

        template<class Cmp, class UpdateCell>
        void sort_with_updates(const Cmp& cmp, UpdateCell update_cell)
        {
            get_order().sort(cmp);
            update_indices(update_cell);
        }

        template<class Cmp, class UpdateCell>
        void sort_with_index_map(const Cmp& cmp, UpdateCell update_cell)
        {
            get_order().sort(cmp);
            auto indices = current_to_sorted_indices(get_order(), size());
            update_indices(indices, update_cell);
        }

        template<class UpdateCell>
        void update_indices(UpdateCell update_cell)
        {
            update_order_indices(get_order(), get_complex(),
                                 [this](typename Order::iterator it) { return project_complex(it); },
                                 [&update_cell](IndexedCell& c, size_t i) { update_cell(c, i); });
        }

        template<class UpdateCell>
        void update_indices(const std::vector<size_t>& indices, UpdateCell update_cell)
        {
            update_order_indices(get_order(), get_complex(),
                                 [this](typename Order::iterator it) { return project_complex(it); },
                                 [&indices,&update_cell](IndexedCell& c, size_t i) { update_cell(c, i, indices); });
        }

    private:
        Container           cells_;
};

}

}

#endif
