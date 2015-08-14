#include <boost/range/adaptors.hpp>
namespace ba = boost::adaptors;

template<class P>
template<class Filtration, class ReportPair>
void
dionysus::StandardReduction<P>::
operator()(const Filtration& filtration, const ReportPair& report_pair)
{
    persistence_.reserve(filtration.size());

    typedef     typename Filtration::Cell                       Cell;
    typedef     ChainEntry<Field, Cell>                         CellChainEntry;
    typedef     ChainEntry<Field, Index>                        ChainEntry;

    unsigned i = 0;
    for(auto& c : filtration)
    {
        //std::cout << "Adding: " << c << " : " << boost::distance(c.boundary(persistence_.field())) << std::endl;
        Index pair = persistence_.add(c.boundary(persistence_.field()) |
                                                 ba::transformed([this,&filtration](const CellChainEntry& e)
                                                 { return ChainEntry(e.element(), filtration.index(e.index())); }));
        if (pair != persistence_.unpaired())
            report_pair(c.dimension(), pair, i);
        ++i;
    }
}
