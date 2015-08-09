#ifndef DIONYSUS_ROW_REDUCTION_H
#define DIONYSUS_ROW_REDUCTION_H

namespace dionysus
{

// Mid-level interface
template<class Field_, typename Index_ = unsigned, class Comparison_ = std::less<Index_>, template<class Persistence> class... Visitors>
class RowReduction
{
    public:
        typedef         Field_                                                  Field;
        typedef         Index_                                                  Index;
        typedef         Comparison_                                             Comparison;

        typedef         ReducedMatrix<Field_,Index_,Comparison_,Visitors...>    Persistence;

    public:
                        RowReduction(const Field& field):
                            persistence_(field)                         {}

                        RowReduction(const Field&                       field,
                                     const Comparison&                  cmp,
                                     const Visitors<Persistence>&...    visitors):
                            persistence_(field, cmp, visitors...)       {}

        template<class Filtration>
        void            operator()(const Filtration& f);

        const Persistence&
                        persistence() const                         { return persistence_; }
        Persistence&    persistence()                               { return persistence_; }

    private:
        Persistence     persistence_;
};

}

#include "row-reduction.hpp"

#endif

