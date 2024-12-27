#ifndef DIONYSUS_MULTI_FILTRATION_H
#define DIONYSUS_MULTI_FILTRATION_H

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

// MultiFiltration stores a filtered cell complex as boost::multi_index_container<...>.
// It allows for bidirectional translation between a cell and its index.
template<class Cell_,
         bool  checked_index = false>
class MultiFiltration
{
    public:
        struct order {};

        typedef             Cell_                                               Cell;

        struct CellWithIndex: Cell
        {
                        CellWithIndex(const Cell& c_, size_t i_):
                            Cell(c_), i(i_)             {}
                        CellWithIndex(Cell&& c_, size_t i_):
                            Cell(c_), i(i_)             {}

            size_t  i;

            friend bool operator<(const CellWithIndex& c, const CellWithIndex& other)
            {
                const Cell& cc = c;
                const Cell& oc = other;
                return cc < oc || (cc == oc && c.i < other.i);
            }
        };
        struct CellWithoutIndex
        {
            Cell&           operator()(CellWithIndex& c) const            { return c; }
            const Cell&     operator()(const CellWithIndex& c) const      { return c; }
        };
        // non-unique to avoid modification conflicts in update_indices()
        using CellLookupIndex = bmi::ordered_non_unique<bmi::identity<CellWithIndex>>;
        using Container = b::multi_index_container<CellWithIndex, bmi::indexed_by<CellLookupIndex,
                                                                                   bmi::random_access<bmi::tag<order>>
                                                                                  >>;
        typedef             typename Container::value_type                      value_type;

        typedef             typename Container::template nth_index<0>::type     Complex;
        typedef             typename Container::template nth_index<1>::type     Order;

        using OrderConstIterator = b::transform_iterator<CellWithoutIndex, typename Order::const_iterator>;
        using OrderIterator = b::transform_iterator<CellWithoutIndex, typename Order::iterator>;


    public:
                            MultiFiltration()                                        = default;
                            MultiFiltration(MultiFiltration&& other)                 = default;
        MultiFiltration&    operator=(MultiFiltration&& other)                       = default;

                            MultiFiltration(const std::initializer_list<Cell>& cells):
                                MultiFiltration(std::begin(cells), std::end(cells))  {}

        template<class Iterator>
                            MultiFiltration(Iterator bg, Iterator end)               { for (auto it = bg; it != end; ++it) cells_.push_back(*it); }

        template<class CellRange>
                            MultiFiltration(const CellRange& cells):
                                MultiFiltration(std::begin(cells), std::end(cells))  {}

        // Lookup
        const Cell&         operator[](size_t i) const                          { return get_order()[i]; }
        size_t              index(const Cell& s, size_t i) const;
        bool                contains(const Cell& s) const                       { return get_complex().contains(s, std::less<Cell>()); }

        void                push_back(const Cell& s)                            { get_order().push_back( CellWithIndex(s, cells_.size()) ); }
        void                push_back(Cell&& s)                                 { get_order().push_back( CellWithIndex(s, cells_.size()) ); }

        void                replace(size_t i, const Cell& s)                    { get_order().replace(get_order().begin() + i, CellWithIndex(s, i)); }

        template<class... Args>
        void                emplace_back(Args&&... args)                        { get_order().emplace_back( CellWithIndex(Cell(std::forward<Args>(args)...), cells_.size()) ); }

        template<class Cmp = std::less<Cell>>
        void                sort(const Cmp& cmp = Cmp())                        { get_order().sort(cmp); update_indices(); }

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

        void                update_indices();

        Container           cells_;
};

}

template<class C, bool checked_index>
size_t
dionysus::MultiFiltration<C,checked_index>::
index(const Cell& s, size_t i) const
{
    auto it = get_complex().upper_bound(CellWithIndex(s,i));
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
dionysus::MultiFiltration<C,checked_index>::
rearrange(const std::vector<size_t>& indices)
{
    std::vector<std::reference_wrapper<const CellWithIndex>> references; references.reserve(indices.size());
    for (size_t i : indices)
    {
        auto& c = get_order()[i];
        references.push_back(std::cref(c));
    }
    get_order().rearrange(references.begin());

    update_indices();
}

template<class C, bool checked_index>
void
dionysus::MultiFiltration<C,checked_index>::
update_indices()
{
    size_t i = 0;
    for(auto it = get_order().begin(); it != get_order().end(); ++it)
    {
        auto cit = project_complex(it);       // complex iterator
        get_complex().modify(cit, [i](CellWithIndex& c) { c.i = i; });
        ++i;
    }
}

#endif
