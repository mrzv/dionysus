#include <iostream>
#include <vector>

#include <dionysus/simplex.h>
#include <dionysus/filtration.h>
#include <dionysus/fields/zp.h>
#include <dionysus/fields/z2.h>
#include <dionysus/ordinary-persistence.h>
#include <dionysus/standard-reduction.h>
#include <dionysus/row-reduction.h>
#include <dionysus/cohomology-persistence.h>
#include <dionysus/zigzag-persistence.h>

namespace d = dionysus;

#include <format.h>

//typedef     d::Z2Field                  K;
typedef     d::ZpField<>                K;
typedef     d::Simplex<>                Simplex;
typedef     d::Filtration<Simplex>      Filtration;
//typedef     d::OrdinaryPersistence<K>   Persistence;
//typedef     d::OrdinaryPersistenceNoNegative<K>   Persistence;
//typedef     d::CohomologyPersistence<K>     Persistence;
typedef     d::ZigzagPersistence<K>       Persistence;

int main()
{
    //K k;
    K k(11);

    Simplex         s {0,1,3};
    for (auto sb : s.boundary())
        fmt::print("{}\n", sb);


    Filtration filtration { Simplex{0}, Simplex{1}, Simplex{2}, Simplex{0,1}, Simplex{0,2}, Simplex{1,2}, Simplex{0,1,2} };

    for (auto& s : filtration)
    {
        fmt::print("{} at {}\n", s, filtration.index(s));
        for (auto sb : s.boundary(k))
            fmt::print("   {} * {} at {}\n", sb.element(), sb.index(), filtration.index(sb.index()));
    }

    //Persistence                             persistence(k);
    //d::StandardReduction<Persistence>       reduce(persistence);
    d::RowReduction<K>                      reduce(k);
    reduce(filtration);
    fmt::print("Reduction finished\n");

#if 0
    unsigned i = 0;
    for (auto& c : persistence.columns())
    {
        std::cout << i << ": ";
        for (auto& ce : c.chain)
        {
            std::cout << " + " << ce.element() << " * " << ce.index();
        }
        std::cout << std::endl;
        ++i;
    }
#endif
}
