template<class F, class I, class C>
template<class ChainRange>
typename dionysus::ZigzagPersistence<F,I,C>::Index
dionysus::ZigzagPersistence<F,I,C>::
add(const ChainRange& chain)
{
    IndexChain cycles;      // chain -> Z*cycles
    Column     z_remainder = Z.reduce(chain, cycles);
    assert(z_remainder.empty());

    IndexChain  boundaries;
    Column      b_remainder = B.reduce(cycles, boundaries);

    // add up columns of D indexed by boundaries
    typedef     typename IndexChain::value_type         IndexChainEntry;
    auto entry_cmp = [this](const IndexChainEntry& e1, const IndexChainEntry& e2)
                     { return this->cmp()(e1.index(), e2.index()); };
    typedef std::set<IndexChainEntry, decltype(entry_cmp)>       IndexChainSet;
    IndexChainSet       boundary(entry_cmp);
    for (auto& x : boundaries)
        Chain<IndexChainSet>::addto(boundary, x.element(), D[x.index()], field(), cmp());
    Index cell = cell_indices++;
    boundary.insert(IndexChainEntry(field().id(), cell));

    if (b_remainder.empty())        // birth
    {
        Index z_col = z_indicies_last++;
        Z.set(z_col, boundary);
        birth_index[z_col] = cell;
        return unpaired;
    }
    else                            // death
    {
        Index b_col = b_indices++;
        Index pair  = std::get<0>(b_remainder.back().index());
        B.set(b_col, std::move(b_remainder));

        IndexChain& column = D[b_col];
        for (auto& x : boundary)
            column.push_back(x);

        return birth_index[pair];
    }
}
