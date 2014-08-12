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

template<class Field_, class Index_ = int, class Comparison_ = std::less<Index_>>
class SparseRowMatrix
{
    public:
        typedef         Field_                      Field;
        typedef         Index_                      Index;
        typedef         Comparison_                 Comparison;

        typedef         typename Field::Element     FieldElement;
        typedef         std::tuple<Index, Index>    IndexPair;                  // (id, pair)

        struct          Entry;                                                  // (FieldElement, (Index, Index))
        typedef         bi::list_base_hook<bi::link_mode<bi::auto_unlink>>      auto_unlink_hook;
        typedef         std::vector<Entry>                                      Column;
        typedef         bi::list<Entry, bi::constant_time_size<false>>          Row;

        typedef         std::vector<ChainEntry<Field, Index>>                   IndexChain;

        typedef         std::unordered_map<Index, Column>                       Columns;
        typedef         std::unordered_map<Index, Row>                          Rows;
        typedef         std::unordered_map<Index, Index>                        LowMap;

    public:
                        SparseRowMatrix(const Field&          field,
                                        const Comparison&     cmp = Comparison()):
                            field_(field), cmp_(cmp)                            {}

        // TODO
        template<class ChainRange>
        Column          reduce(const ChainRange& chain, IndexChain& trail);

        template<class ChainRange>
        void            set(Index i, const ChainRange& chain);

        void            set(Index i, Column&& chain);

        // TODO: remove

        // accessors
        Row&            row(Index r)                                            { return rows_[r]; }
        Column&         col(Index c)                                            { return columns_[c]; }
        Index           low(Index r) const                                      { return lows_[r]; }            // column that has this low

        const Field&    field() const                                           { return field_; }
        void            reserve(size_t)                                         {}                              // here for compatibility only
        const Comparison&   cmp() const                                         { return cmp_; }

    private:
        Field       field_;
        Comparison  cmp_;

        Columns     columns_;
        Rows        rows_;
        LowMap      lows_;
};

template<class F, class I, class C>
struct SparseRowMatrix<F,I,C>::Entry:
    public ChainEntry<F, IndexPair, auto_unlink_hook>
{
    typedef             ChainEntry<Field, IndexPair, auto_unlink_hook>      Parent;

                        Entry(FieldElement e, const IndexPair& ip):
                            Parent(e,ip)                                    {}

                        Entry(FieldElement e, const Index& r, const Index& c):
                            Parent(e,IndexPair(r,c))                        {}

                        Entry(const Entry& other)   = default;
                        Entry(Entry&& other)        = default;
    Entry&              operator=(Entry&& other)    = default;

    void                unlink()                                            { auto_unlink_hook::unlink(); }
    bool                is_linked()  const                                  { return auto_unlink_hook::is_linked();  }
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
