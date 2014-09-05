#include <iostream>
#include <fstream>

#include <dionysus/rips.h>
#include <dionysus/filtration.h>
#include <dionysus/fields/zp.h>
#include <dionysus/fields/z2.h>
#include <dionysus/ordinary-persistence.h>
#include <dionysus/standard-reduction.h>

namespace d = dionysus;

// Trivial example of size() points on a line with integer coordinates
struct Distances
{
    typedef         int             IndexType;
    typedef         double          DistanceType;

    DistanceType    operator()(IndexType a, IndexType b) const      { return std::abs(a - b); }

    size_t          size() const                                    { return 2000; }
    IndexType       begin() const                                   { return 0; }
    IndexType       end() const                                     { return size(); }
};

//typedef         Rips<ExplicitDistances<Distances> >                   Generator;
typedef         d::Rips<Distances>                                      Generator;
typedef         Generator::Simplex                                      Simplex;
typedef         d::Filtration<Simplex>                                  Filtration;
//typedef       d::Z2Field                                                K;
typedef         d::ZpField<>                                            K;
typedef         d::OrdinaryPersistence<K>                               Persistence;


int main(int argc, char* argv[])
{
    Distances distances;

    // Storing ExplicitDistances speeds up the computation (at the price of memory)
    //ExplicitDistances<Distances> explicit_distances(distances);

    Generator               rips(distances);
    Generator::Evaluator    size(distances);
    Filtration              filtration;

    // Generate 2-skeleton of the Rips complex for epsilon = 50
    rips.generate(2, 10,
                  [&filtration](Simplex&& s) { filtration.push_back(s); });
    std::cout << "Generated complex of size: " << filtration.size() << std::endl;

    // Generate filtration with respect to distance and compute its persistence
    filtration.sort(Generator::Comparison(distances));

    //K k;
    K k(11);
    Persistence                             persistence(k);
    d::StandardReduction<Persistence>       reduce(persistence);
    reduce(filtration);
    std::cout << "Reduction finished" << std::endl;
}
