template<class F, class I, class C>
template<class ChainRange>
typename dionysus::SparseRowMatrix<F,I,C>::Column
dionysus::SparseRowMatrix<F,I,C>::
reduce(const ChainRange& chain_, IndexChain& trail)
{
    auto    row_cmp = [this](const Entry& e1, const Entry& e2)
                      { return this->cmp_(std::get<0>(e1.index()), std::get<0>(e2.index())); };

#define __DIONYSUS_USE_VECTOR_CHAINS    0

#if !(__DIONYSUS_USE_VECTOR_CHAINS)
    std::set<Entry,decltype(row_cmp)>   chain(row_cmp);
    for (auto x : chain_)
        chain.insert(Entry(x.element(), IndexPair(x.index(), 0)));
#else
    Column chain;
    for (auto x : chain_)
        chain.emplace_back(x.element(), IndexPair(x.index(), 0));
    std::sort(chain.begin(), chain.end(), row_cmp);
#endif

    typedef   Reduction<IndexPair>              ReductionIP;

    auto      chains   = [this](const IndexPair& rc) -> const Column&    { return this->col(std::get<1>(rc)); };
    auto      lows     = [this](const IndexPair& rc) -> IndexPair
                         {
                             Index r  = std::get<0>(rc);
                             Index c  = std::get<1>(rc);
                             auto  it = this->lows_.find(r);
                             if (it == this->lows_.end())
                                 return ReductionIP::unpaired;
                             else
                                 return IndexPair(r, it->second);
                         };

    auto      addto    = [&trail](FieldElement m, const IndexPair& rc)  { trail.emplace_back(m, std::get<1>(rc)); };

    IndexPair pair     = ReductionIP::reduce(chain,
                                             chains, lows,
                                             field_, addto, row_cmp);

#if !(__DIONYSUS_USE_VECTOR_CHAINS)
    return Column(std::begin(chain), std::end(chain));
#else
    return chain;
#endif
}

template<class F, class I, class C>
template<class ChainRange>
void
dionysus::SparseRowMatrix<F,I,C>::
set(Index col, const ChainRange& chain)
{
    Column& column = columns_.emplace(col, Column()).first->second;
    column.reserve(chain.size());

    for (auto& x : chain)
    {
        Index r = x.index();
        column.push_back(Entry(x.element(), IndexPair(r, col)));
        row(r).push_back(column.back());
    }

    Index r = std::get<0>(column.back().index());
    lows_[r] = col;
}

template<class F, class I, class C>
void
dionysus::SparseRowMatrix<F,I,C>::
set(Index col, Column&& chain)
{
    Column& column = columns_.emplace(col, std::move(chain)).first->second;

    for (auto& x : column)
    {
        std::get<1>(x.index()) = col;
        Index r = std::get<0>(x.index());
        row(r).push_back(x);
    }

    Index r = std::get<0>(column.back().index());
    lows_[r] = col;
}
