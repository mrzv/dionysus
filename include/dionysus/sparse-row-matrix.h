#ifndef DIONYSUS_SPARSE_ROW_MATRIX_H
#define DIONYSUS_SPARSE_ROW_MATRIX_H

#include <vector>
#include <list>
#include <unordered_map>

#include <boost/intrusive/list.hpp>

#include "chain.h"
#include "reduction.h"

namespace dionysus
{

namespace bi = boost::intrusive;

namespace detail
{
    typedef         bi::list_base_hook<bi::link_mode<bi::auto_unlink>>      auto_unlink_hook;

    template<class F, class I>
    struct SparseRowMatrixEntry:
        public ChainEntry<F, std::tuple<I,I>, auto_unlink_hook>
    {
        typedef             I                                                   Index;
        typedef             typename F::Element                                 FieldElement;
        typedef             std::tuple<Index, Index>                            IndexPair;                  // (id, pair)
        typedef             ChainEntry<F, IndexPair, auto_unlink_hook>          Parent;
        typedef             SparseRowMatrixEntry                                Entry;

                            SparseRowMatrixEntry(FieldElement e, const IndexPair& ip):
                                Parent(e,ip)                                    {}

                            SparseRowMatrixEntry(FieldElement e, const Index& r, const Index& c):
                                Parent(e,IndexPair(r,c))                        {}

                            SparseRowMatrixEntry(const Entry& other)   = default;
                            SparseRowMatrixEntry(Entry&& other)        = default;
        Entry&              operator=(Entry&& other)    = default;

        void                unlink()                                            { auto_unlink_hook::unlink(); }
        bool                is_linked()  const                                  { return auto_unlink_hook::is_linked();  }
    };
}

template<class Field_, class Index_ = int, class Comparison_ = std::less<Index_>,
         template<class E, class... Args> class Column_ = std::vector>
class SparseRowMatrix
{
    public:
        typedef         Field_                                                  Field;
        typedef         Index_                                                  Index;
        typedef         Comparison_                                             Comparison;

        typedef         typename Field::Element                                 FieldElement;

        typedef         detail::SparseRowMatrixEntry<Field,Index>               Entry;
        typedef         Column_<Entry>                                          Column;
        typedef         typename Entry::IndexPair                               IndexPair;
        typedef         bi::list<Entry, bi::constant_time_size<false>>          Row;

        typedef         std::vector<ChainEntry<Field, Index>>                   IndexChain;

        typedef         std::unordered_map<Index, Column>                       Columns;
        typedef         std::unordered_map<Index, Row>                          Rows;
        typedef         std::unordered_map<Index, Index>                        LowMap;

    public:
                        SparseRowMatrix(const Field&        field,
                                        const Comparison&   cmp = Comparison()):
                            field_(field), cmp_(cmp)                            {}

        template<class ChainRange>
        Column          reduce(const ChainRange& chain, IndexChain& trail);

        Index           set(Index i, Column&& chain);       // returns previous column with this low
        void            fix(Index c, Column& column);
        void            fix(Index c)                                            { fix(c, col(c)); }

        void            prepend_row(Index r, FieldElement m, const Row& chain); // could be horribly inefficient if Column is chosen poorly

        void            drop_row(Index r)                                       { rows_.erase(r); if (is_low(r)) lows_.erase(r); }
        void            drop_col(Index c)                                       { columns_.erase(c); }

        // accessors
        Row&            row(Index r)                                            { return rows_[r]; }
        Column&         col(Index c)                                            { return columns_[c]; }
        const Column&   col(Index c) const                                      { return columns_.find(c)->second; }
        Index           low(Index r) const                                      { return lows_.find(r)->second; }
        bool            is_low(Index r) const                                   { return lows_.find(r) != lows_.end(); }
        void            update_low(Index c)                                     { lows_[std::get<0>(col(c).back().index())] = c; }

        const Field&        field() const                                       { return field_; }
        void                reserve(size_t)                                     {}                              // here for compatibility only
        const Comparison&   cmp() const                                         { return cmp_; }

    private:
        Field           field_;
        Comparison      cmp_;

        Columns         columns_;
        Rows            rows_;
        LowMap          lows_;          // column that has this low
};


namespace detail
{

template<class Index>
struct Unpaired<std::tuple<Index,Index>>
{
    static
    constexpr std::tuple<Index,Index>
    value()
    { return std::make_tuple(std::numeric_limits<Index>::max(),
                             std::numeric_limits<Index>::max()); }
};

}

}

#include "sparse-row-matrix.hpp"

#endif
