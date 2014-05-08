#include <boost/range/adaptors.hpp>
namespace ba = boost::adaptors;

template<class P>
template<class Filtration>
void
dionysus::StandardReduction<P>::
operator()(const Filtration& filtration)
{
    persistence_.reserve(filtration.size());

    typedef     ChainEntry<Field, typename Filtration::Cell>      CellChainEntry;
    typedef     ChainEntry<Field, typename Persistence::Index>    ChainEntry;

    for(auto& c : filtration)
    {
        std::cout << "Adding: " << c << std::endl;
        persistence_.add(c.boundary(persistence_.field()) |
                         ba::transformed([this,&filtration](const CellChainEntry& e)
                         { return ChainEntry(e.element(), filtration.index(e.index())); }));
    }
}
