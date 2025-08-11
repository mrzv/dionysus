#ifndef DIONYSUS_LINKED_MULTI_FILTRATION_H
#define DIONYSUS_LINKED_MULTI_FILTRATION_H

#include <vector>
#include <sstream>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/random_access_index.hpp>

#include <boost/iterator/transform_iterator.hpp>

namespace b   = boost;
namespace bmi = boost::multi_index;

namespace dionysus
{

// LinkedMultiFiltration stores for each cell a linked index. The linked index is checked first by the index(c,i) search.
// TODO: there is a lot of code duplication here with MultiFiltration.
template<class Cell_,
         bool  checked_index = false>
class LinkedMultiFiltration
{
    public:
        struct order {};

        typedef             Cell_                                               Cell;

        struct LinkedCellWithIndex: Cell
        {
                        LinkedCellWithIndex(const Cell& c_, size_t i_, size_t l_):
                            Cell(c_), i(i_), linked(l_)             {}
                        LinkedCellWithIndex(Cell&& c_, size_t i_, size_t l_):
                            Cell(c_), i(i_), linked(l_)             {}

            size_t  i;
            size_t  linked;

            friend bool operator<(const LinkedCellWithIndex& c, const LinkedCellWithIndex& other)
            {
                const Cell& cc = c;
                const Cell& oc = other;
                return cc < oc || (cc == oc && c.i < other.i);
            }
        };
        struct CellWithoutIndex
        {
            Cell&           operator()(LinkedCellWithIndex& c) const            { return c; }
            const Cell&     operator()(const LinkedCellWithIndex& c) const      { return c; }
        };
        // non-unique to avoid modification conflicts in update_indices()
        using CellLookupIndex = bmi::ordered_non_unique<bmi::identity<LinkedCellWithIndex>>;
        using Container = b::multi_index_container<LinkedCellWithIndex,
                                                   bmi::indexed_by<CellLookupIndex,
                                                                   bmi::random_access<bmi::tag<order>>
                                                                  >>;
        typedef             typename Container::value_type                      value_type;

        typedef             typename Container::template nth_index<0>::type     Complex;
        typedef             typename Container::template nth_index<1>::type     Order;

        using OrderConstIterator = b::transform_iterator<CellWithoutIndex, typename Order::const_iterator>;
        using OrderIterator = b::transform_iterator<CellWithoutIndex, typename Order::iterator>;


    public:
                            LinkedMultiFiltration()                                        = default;
                            LinkedMultiFiltration(LinkedMultiFiltration&& other)           = default;
        LinkedMultiFiltration&
                            operator=(LinkedMultiFiltration&& other)                       = default;

                            LinkedMultiFiltration(const std::initializer_list<Cell>& cells):
                                LinkedMultiFiltration(std::begin(cells), std::end(cells))  {}

        template<class Iterator>
                            LinkedMultiFiltration(Iterator bg, Iterator end)               { for (auto it = bg; it != end; ++it) cells_.push_back(*it); }

        template<class CellRange>
                            LinkedMultiFiltration(const CellRange& cells):
                                LinkedMultiFiltration(std::begin(cells), std::end(cells))  {}

        // Lookup
        const Cell&         operator[](size_t i) const                          { return get_order()[i]; }
        size_t              index(const Cell& s, size_t i) const;
        bool                contains(const Cell& s) const                       { return get_complex().contains(s, std::less<Cell>()); }

        void                push_back(const Cell& s, size_t l)                  { get_order().push_back( LinkedCellWithIndex(s, cells_.size(), l) ); }
        void                push_back(Cell&& s, size_t l)                       { get_order().push_back( LinkedCellWithIndex(s, cells_.size(), l) ); }

        void                replace(size_t i, const Cell& s, size_t l)          { get_order().replace(get_order().begin() + i, LinkedCellWithIndex(s, i, l)); }

        template<class... Args>
        void                emplace_back(size_t l, Args&&... args)              { get_order().emplace_back( LinkedCellWithIndex(Cell(std::forward<Args>(args)...), cells_.size(), l) ); }

        template<class Cmp = std::less<Cell>>
        void                sort(const Cmp& cmp = Cmp());

        void                rearrange(const std::vector<size_t>& indices);

        OrderConstIterator  begin() const                                       { return OrderConstIterator(get_order().begin()); }
        OrderConstIterator  end() const                                         { return OrderConstIterator(get_order().end()); }
        OrderIterator       begin()                                             { return OrderIterator(get_order().begin()); }
        OrderIterator       end()                                               { return OrderIterator(get_order().end()); }
        size_t              size() const                                        { return cells_.size(); }
        void                clear()                                             { return Container().swap(cells_); }

        Cell&               back()                                              { return const_cast<Cell&>(get_order().back()); }
        const Cell&         back() const                                        { return get_order().back(); }

    private:
        const Complex&      get_complex() const                                 { return cells_.template get<0>(); }
        Complex&            get_complex()                                       { return cells_.template get<0>(); }
        const Order&        get_order() const                                   { return cells_.template get<order>(); }
        Order&              get_order()                                         { return cells_.template get<order>(); }

        template<class Iterator>
        typename Order::const_iterator
                            project_order(Iterator it) const                    { return cells_.template project<order>(it); }
        template<class Iterator>
        typename Complex::iterator
                            project_complex(Iterator it)                        { return cells_.template project<0>(it); }

        void                update_indices(const std::vector<size_t>& indices);

        Container           cells_;
};

}

template<class C, bool checked_index>
size_t
dionysus::LinkedMultiFiltration<C,checked_index>::
index(const Cell& s, size_t i) const
{
    size_t l = get_order()[i].linked;
    auto& cl = (*this)[l];
    if (cl == s)
        return l;

    // linked = 0 because linked is not used in complex order
    auto it = get_complex().upper_bound(LinkedCellWithIndex(s,i,0));
    --it;

    if (checked_index && *it != s)
    {
        std::ostringstream oss;
        oss << "Trying to access non-existent cell: " << s << ' ' << i;
        throw std::runtime_error(oss.str());
    }

    return project_order(it) - get_order().begin();
}

template<class C, bool checked_index>
void
dionysus::LinkedMultiFiltration<C,checked_index>::
rearrange(const std::vector<size_t>& indices)
{
    std::vector<std::reference_wrapper<const LinkedCellWithIndex>> references; references.reserve(indices.size());
    std::vector<size_t> reverse_indices(indices.size());
    size_t i = 0;
    for (size_t idx : indices)
    {
        reverse_indices[idx] = i++;
        auto& c = get_order()[idx];
        references.push_back(std::cref(c));
    }
    get_order().rearrange(references.begin());

    update_indices(reverse_indices);
}

template<class C, bool checked_index>
template<class Cmp>
void
dionysus::LinkedMultiFiltration<C,checked_index>::
sort(const Cmp& cmp)
{
    get_order().sort(cmp);

    // record where the old cells landed
    size_t i = 0;
    std::vector<size_t> indices(size());
    for (auto& c : get_order())
        indices[c.i] = i++;

    update_indices(indices);
}

template<class C, bool checked_index>
void
dionysus::LinkedMultiFiltration<C,checked_index>::
update_indices(const std::vector<size_t>& indices)
{
    size_t i = 0;
    for(auto it = get_order().begin(); it != get_order().end(); ++it)
    {
        auto cit = project_complex(it);       // complex iterator
        get_complex().modify(cit, [i,&indices](LinkedCellWithIndex& c) { c.i = i; c.linked = indices[c.linked]; });
        ++i;
    }
}

#endif
