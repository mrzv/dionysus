#ifndef DIONYSUS_CLEARING_REDUCTION_H
#define DIONYSUS_CLEARING_REDUCTION_H

namespace dionysus
{

// Mid-level interface
template<class Persistence_>
class ClearingReduction
{
    public:
        typedef     Persistence_                Persistence;

    public:
                    ClearingReduction(Persistence& persistence):
                        persistence_(persistence)               {}

        template<class Filtration>
        void        operator()(const Filtration& f);

    private:
        Persistence&  persistence_;
};

}

#include "clearing-reduction.hpp"

#endif

