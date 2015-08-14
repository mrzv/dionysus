#include <dionysus/relative-homology-zigzag.h>
#include <dionysus/fields/zp.h>
#include <dionysus/fields/z2.h>

#include <dionysus/simplex.h>
#include <dionysus/filtration.h>
namespace d = dionysus;

#include <boost/range/adaptors.hpp>
namespace ba = boost::adaptors;

#include <format.h>

typedef     d::Z2Field                      K;
//typedef     d::ZpField<>                    K;
typedef     d::Simplex<>                    Simplex;
typedef     d::RelativeHomologyZigzag<K>    Persistence;

typedef     d::Filtration<Simplex>          Filtration;
typedef     typename Filtration::Cell       Cell;

typedef     typename Persistence::Index     Index;
typedef     d::ChainEntry<K, Cell>          CellChainEntry;
typedef     d::ChainEntry<K, Index>         ChainEntry;


int main()
{
    K k;
    //K k(11);
    Filtration      filtration { Simplex{0}, Simplex{1}, Simplex{2}, Simplex{0,1}, Simplex{0,2}, Simplex{1,2}, Simplex{0,1,2} };
    Persistence     persistence(k);

    unsigned op = 0;
    for(auto& c : filtration)
    {
        fmt::print("[{}] Adding: {} : {}\n", op++, c, boost::distance(c.boundary(persistence.field())));
        persistence.add_both(c.boundary(persistence.field()) |
                                        ba::transformed([&filtration](const CellChainEntry& e)
                                        {
                                            Index idx = filtration.index(e.index());
                                            return ChainEntry(e.element(), idx);
                                        }));
    }

    for (int i = 6; i >= 0; --i)
    {
        fmt::print("[{}] Removing: {} = {}\n", op++, i, filtration[i]);
        Index pair = persistence.remove(i);
        if (pair == Persistence::unpaired())
            fmt::print("Birth\n");
        else
            fmt::print("Death: {}\n", pair);
    }
}

