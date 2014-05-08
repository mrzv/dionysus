#include <iostream>
#include <vector>

#include <dionysus/simplex.h>
#include <dionysus/filtration.h>
#include <dionysus/fields/zp.h>
#include <dionysus/fields/z2.h>
#include <dionysus/ordinary-persistence.h>
#include <dionysus/standard-reduction.h>

namespace d = dionysus;

//typedef     d::Z2Field                  K;
typedef     d::ZpField<>                K;
typedef     d::Simplex<>                Simplex;
typedef     d::Filtration<Simplex>      Filtration;
//typedef     d::OrdinaryPersistence<K>   Persistence;
typedef     d::OrdinaryPersistenceNoNegative<K>   Persistence;

int main()
{
    //K k;
    K k(11);

    Simplex         s {0,1,3};
    for (auto sb : s.boundary())
        std::cout << sb << std::endl;


    Filtration filtration { Simplex{0}, Simplex{1}, Simplex{2}, Simplex{0,1}, Simplex{0,2}, Simplex{1,2}, Simplex{0,1,2} };

    for (auto& s : filtration)
    {
        std::cout << s << " at " << filtration.index(s) << std::endl;
        for (auto sb : s.boundary(k))
            std::cout << "   " << sb.element() << " * " << sb.index() << " at " << filtration.index(sb.index()) << std::endl;
    }

    Persistence                             persistence(k);
    d::StandardReduction<Persistence>       reduce(persistence);
    reduce(filtration);
    std::cout << "Reduction finished" << std::endl;

    unsigned i = 0;
    for (auto& c : persistence.columns())
    {
        std::cout << i << ": ";
        for (auto& ce : c)
        {
            std::cout << " + " << ce.element() << " * " << ce.index();
        }
        std::cout << std::endl;
        ++i;
    }
}
