#ifndef DIONYSUS_REDUCTION_H
#define DIONYSUS_REDUCTION_H

#include <vector>
#include "chain.h"

namespace dionysus
{

template<class Index_>
struct Reduction
{
    typedef     Index_              Index;

    template<class Field>
    using AddtoVisitor = std::function<void(typename Field::Element, Index)>;

    static const Index unpaired;

    template<class Chain1,
             class ChainsLookup,
             class LowLookup,
             class Field,
             class Comparison = std::less<Index>>
    static
    Index reduce(Chain1&                     c,
                 const ChainsLookup&         chains,
                 const LowLookup&            lows,
                 const Field&                field,
                 const AddtoVisitor<Field>&  visitor = [](typename Field::Element, Index) {},
                 const Comparison&           cmp     = Comparison())
    {
        typedef     typename Field::Element         FieldElement;

        while (!c.empty())
        {
            //auto&  low = c.back();
            auto&  low = *(--c.end());
            Index  l   = low.index();
            Index  cl  = lows(l);
            if (cl == unpaired)
                return l;
            else
            {
                // Reduce further
                auto&           co     = chains(cl);
                auto&           co_low = co.back();
                FieldElement    m      = field.neg(field.div(low.element(), co_low.element()));
                // c += m*co
                Chain<Chain1>::addto(c, m, co, field, cmp);
                visitor(m, cl);
            }
        }
        return unpaired;
    }

    template<class Chain1,
             class Chain2,
             class Field,
             class Comparison = std::less<Index>>
    static
    Index reduce(Chain1&                     c,
                 const std::vector<Chain2>&  chains,
                 const std::vector<Index>&   lows,
                 const Field&                field,
                 const AddtoVisitor<Field>&  visitor = [](typename Field::Element, Index) {},
                 const Comparison&           cmp     = Comparison())
    {
        return reduce(c,
                      [&chains](Index i)    { return chains[i]; },
                      [&lows](Index i)      { return lows[i]; },
                      field, visitor, cmp);
    }
};


template<class Index>
const Index
Reduction<Index>::unpaired = std::numeric_limits<Index>::max();

}

#endif
