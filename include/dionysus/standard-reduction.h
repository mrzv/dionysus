#ifndef DIONYSUS_STANDARD_REDUCTION_H
#define DIONYSUS_STANDARD_REDUCTION_H

namespace dionysus
{

// Mid-level interface
template<class Persistence_>
class StandardReduction
{
    public:
        typedef     Persistence_                                Persistence;
        typedef     typename Persistence::Field                 Field;

    public:
                    StandardReduction(Persistence& persistence):
                        persistence_(persistence)           {}

        template<class Filtration>
        void        operator()(const Filtration& f);


    private:
        Persistence&  persistence_;
};

}

#include "standard-reduction.hpp"

#endif
