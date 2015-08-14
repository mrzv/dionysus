#include <iostream>
#include <vector>

#include <boost/range/adaptors.hpp>
namespace ba = boost::adaptors;

#include <dionysus/simplex.h>
#include <dionysus/filtration.h>
#include <dionysus/fields/zp.h>
#include <dionysus/fields/z2.h>
#include <dionysus/zigzag-persistence.h>

namespace d = dionysus;

#include <format.h>

//typedef     d::Z2Field                      K;
typedef     d::ZpField<>                    K;
typedef     d::Simplex<>                    Simplex;
typedef     d::Filtration<Simplex>          Filtration;
//typedef     d::ZigzagFiltration<Simplex>    Filtration;
typedef     d::ZigzagPersistence<K>         Persistence;

typedef     typename Persistence::Index                     Index;
typedef     typename Filtration::Cell                       Cell;
typedef     d::ChainEntry<K, Cell>                          CellChainEntry;
typedef     d::ChainEntry<K, Index>                         ChainEntry;


int main()
{
    K k(11);
    Filtration      filtration { Simplex{0}, Simplex{1}, Simplex{2}, Simplex{0,1}, Simplex{0,2}, Simplex{1,2}, Simplex{0,1,2} };
    Persistence     persistence(k);

    unsigned op = 0;
    for(auto& c : filtration)
    {
        fmt::print("[{}] Adding: {} : {}\n", op++, c, boost::distance(c.boundary(persistence.field())));
        Index pair = persistence.add(c.boundary(persistence.field()) |
                                                ba::transformed([&filtration](const CellChainEntry& e)
                                                { return ChainEntry(e.element(), filtration.index(e.index())); }));
        //if (pair != persistence.unpaired())
        //    std::cout << "[" << pair << " - " << i << "]" << std::endl;
        //++i;
    }

    for (int i = 6; i >= 0; --i)
    {
        fmt::print("[{}] Removing: {}\n", op++, i);
        Index pair = persistence.remove(i);
        if (pair == Persistence::unpaired())
            fmt::print("Birth\n");
        else
            fmt::print("Death: {}\n", pair);
    }
}
