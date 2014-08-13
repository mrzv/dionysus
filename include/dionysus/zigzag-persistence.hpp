template<class F, class I, class C>
template<class ChainRange>
typename dionysus::ZigzagPersistence<F,I,C>::Index
dionysus::ZigzagPersistence<F,I,C>::
add(const ChainRange& chain_)
{
    Index op = operations++;

    IndexChain cycles;      // chain_ -> Z*cycles
    Column     z_remainder = Z.reduce(chain_, cycles);
    assert(z_remainder.empty());

    IndexChain  boundaries;
    DequeColumn b_remainder = B.reduce(cycles, boundaries);

    // add up columns of D indexed by boundaries
    typedef     typename Column::value_type             Entry;
    auto        row_cmp = [this](const Entry& e1, const Entry& e2)
                          { return this->cmp()(std::get<0>(e1.index()), std::get<0>(e2.index())); };
    Column      chain;
    for (auto& x : boundaries)
        Chain<Column>::addto(chain, x.element(), D.col(x.index()), field(), row_cmp);
    chain.push_back(Entry(field().id(), IndexPair(cell_indices++,0)));

    if (b_remainder.empty())        // birth
    {
        Index z_col = z_indicies_last++;
        Z.set(z_col, std::move(chain));
        birth_index[z_col] = op;
        return unpaired;
    }
    else                            // death
    {
        Index b_col = b_indices++;
        Index pair  = std::get<0>(b_remainder.back().index());
        B.set(b_col, std::move(b_remainder));
        D.set(b_col, std::move(chain));
        return birth_index[pair];
    }
}
