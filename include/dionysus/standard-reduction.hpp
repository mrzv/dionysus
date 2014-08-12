#include <boost/range/adaptors.hpp>
namespace ba = boost::adaptors;

template<class P>
template<class Filtration>
void
dionysus::StandardReduction<P>::
operator()(const Filtration& filtration)
{
    persistence_.reserve(filtration.size());

    typedef     typename Persistence::Index                     Index;
    typedef     typename Filtration::Cell                       Cell;
    typedef     ChainEntry<Field, Cell>                         CellChainEntry;
    typedef     ChainEntry<Field, Index>                        ChainEntry;

    unsigned i = 0;
    for(auto& c : filtration)
    {
        std::cout << "Adding: " << c << " : " << boost::distance(c.boundary(persistence_.field())) << std::endl;
        Index pair = persistence_.add(c.boundary(persistence_.field()) |
                                                 ba::transformed([this,&filtration](const CellChainEntry& e)
                                                 { return ChainEntry(e.element(), filtration.index(e.index())); }));
        if (pair != persistence_.unpaired)
            std::cout << "[" << pair << " - " << i << "]" << std::endl;
        ++i;
    }
}
