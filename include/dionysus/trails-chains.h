#ifndef DIONYSUS_TRAILS_CHAINS_H
#define DIONYSUS_TRAILS_CHAINS_H

#include "ordinary-persistence.h"

namespace dionysus
{

template<class Field, class Index>
struct MatrixV
{
    template<class Self>
    struct Visitor: public EmptyVisitor<Field, Index, Self>
    {
        void        resized(Self* r, Index sz)                              { v_.resize(sz); }

        template<class Chain>
        void        chain_initialized(Self* r, Index i, Chain& c)           { v_[i].emplace_back(r->field().id(), i); }

        void        addto(Self* r, Index i, typename Field::Element m, Index o)
        {
            Chain<typename Self::Chain>::addto(v_[i], m, v_[o], r->field(), entry_cmp);
        }

        using Entry = typename Self::Entry;
        static bool entry_cmp(const Entry& e1, const Entry& e2)     { return e1.index() < e2.index(); }

        using Comparison = std::less<Index>;
        Comparison  cmp_;

        using Chains = typename Self::Chains;
        Chains      v_;
    };
};

template<class Field, class Index>
struct MatrixU
{
    template<class Self>
    struct Visitor: public EmptyVisitor<Field, Index, Self>
    {
        void        resized(Self* r, Index sz)                              { u_.resize(sz); }

        template<class Chain>
        void        chain_initialized(Self* r, Index i, Chain& c)           { u_[i].emplace_back(r->field().id(), i); }

        void        addto(Self* r, Index i, typename Field::Element m, Index o)
        {
            // Since i is the column being reduced, row i of U is still the
            // identity row. Updating U_o -= m U_i is therefore one append.
            u_[o].emplace_back(r->field().neg(m), i);
        }

        using Chains = typename Self::Chains;
        Chains      u_;
    };
};

template<class    Field,
         typename Index = unsigned,
         class    Comparison = std::less<Index>,
         template<class Self> class... Visitors>
using OrdinaryPersistenceWithV = OrdinaryPersistence<Field, Index, Comparison, MatrixV<Field,Index>::template Visitor, Visitors...>;

template<class    Field,
         typename Index = unsigned,
         class    Comparison = std::less<Index>,
         template<class Self> class... Visitors>
using OrdinaryPersistenceWithU = OrdinaryPersistence<Field, Index, Comparison, MatrixU<Field,Index>::template Visitor, Visitors...>;

template<class    Field,
         typename Index = unsigned,
         class    Comparison = std::less<Index>,
         template<class Self> class... Visitors>
using OrdinaryPersistenceNoNegativeWithV = OrdinaryPersistenceNoNegative<Field, Index, Comparison, MatrixV<Field,Index>::template Visitor, Visitors...>;

template<class    Field,
         typename Index = unsigned,
         class    Comparison = std::less<Index>,
         template<class Self> class... Visitors>
using OrdinaryPersistenceNoNegativeWithU = OrdinaryPersistenceNoNegative<Field, Index, Comparison, MatrixU<Field,Index>::template Visitor, Visitors...>;

}

#endif
