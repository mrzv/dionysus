#ifndef DIONYSUS_FILTRATION_H
#define DIONYSUS_FILTRATION_H

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/random_access_index.hpp>

namespace b   = boost;
namespace bmi = boost::multi_index;

namespace dionysus
{

template<class Cell_,
         class CellLookupIndex_ = bmi::hashed_unique<bmi::identity<Cell_>>>
class Filtration
{
    public:
        struct order {};

        typedef             Cell_                                               Cell;
        typedef             CellLookupIndex_                                    CellLookupIndex;

        typedef             b::multi_index_container<Cell, bmi::indexed_by<CellLookupIndex,
                                                                           bmi::random_access<bmi::tag<order>>
                                                                          >>    Container;
        typedef             typename Container::value_type                      value_type;

        typedef             typename Container::template nth_index<0>::type     Complex;
        typedef             typename Container::template nth_index<1>::type     Order;
        typedef             typename Order::const_iterator                      OrderConstIterator;


    public:
                            Filtration(const std::initializer_list<Cell>& cells):
                                Filtration(std::begin(cells), std::end(cells))  {}

        template<class Iterator, class Cmp = std::less<Cell>>
                            Filtration(Iterator bg, Iterator end,
                                       const Cmp& cmp = Cmp()):
                                cells_(bg, end)                                 { sort(cmp); }

        template<class CellRange, class Cmp = std::less<Cell>>
                            Filtration(const CellRange& cells, const Cmp& cmp = Cmp()):
                                Filtration(std::begin(cells), std::end(cells), cmp) {}

        // Lookup
        const Cell&         operator[](size_t i) const                          { return cells_.template get<order>()[i]; }
        OrderConstIterator  iterator(const Cell& s) const                       { return bmi::project<order>(cells_, cells_.find(s)); }
        size_t              index(const Cell& s) const                          { return iterator(s) - begin(); }

        void                push_back(const Cell& s)                            { cells_.template get<order>().push_back(s); }

        template<class Cmp = std::less<Cell>>
        void                sort(const Cmp& cmp = Cmp())                        { cells_.template get<order>().sort(cmp); }

        OrderConstIterator  begin() const                                       { return cells_.template get<order>().begin(); }
        OrderConstIterator  end() const                                         { return cells_.template get<order>().end(); }
        size_t              size() const                                        { return cells_.size(); }


    private:
        Container           cells_;
};

}

#endif
