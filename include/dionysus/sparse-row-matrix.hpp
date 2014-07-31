template<class F, class I, class C>
template<class ChainRange, class AddTo>
typename dionysus::SparseRowMatrix<F,I,C>::IndexPair
dionysus::SparseRowMatrix<F,I,C>::
add(const ChainRange& chain_, const AddTo& addto)
{
    Index   col     = column_last_++;
    Column& column  = columns_.emplace(col, Column()).first->second;

    auto row_cmp    = [this](const Entry& e1, const Entry& e2)
                      { return this->cmp_(std::get<0>(e1.index()), std::get<0>(e2.index())); };
    std::set<Entry,decltype(row_cmp)>     chain(row_cmp);
    for (auto x : chain_)
        chain.insert(Entry(x.element(), IndexPair(x.index(), col)));

    typedef   Reduction<IndexPair>              Reduction;

    auto      chains   = [this](const IndexPair& rc) -> const Column&    { return this->col(std::get<1>(rc)); };
    auto      lows     = [this](const IndexPair& rc) -> IndexPair
                         {
                             Index r  = std::get<0>(rc);
                             Index c  = std::get<1>(rc);
                             auto  it = this->lows_.find(r);
                             if (it == this->lows_.end())
                                 return Reduction::unpaired;
                             else
                                 return IndexPair(r, it->second);
                         };
    IndexPair pair     = Reduction::reduce(chain,
                                           chains, lows,
                                           field_, addto, row_cmp);

    if (pair != Reduction::unpaired)
        lows_[std::get<0>(pair)] = col;

    column.reserve(chain.size());
    for (auto& x : chain)
    {
        column.push_back(Entry(x.element(), IndexPair(std::get<0>(x.index()), col)));
        Index r = std::get<0>(x.index());
        row(r).push_back(column.back());
    }

    return pair;
}
