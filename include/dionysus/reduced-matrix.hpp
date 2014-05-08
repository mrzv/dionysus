template<class F, typename I, class C, template<class Self> class... V>
template<class ChainRange>
typename dionysus::ReducedMatrix<F,I,C,V...>::Index
dionysus::ReducedMatrix<F,I,C,V...>::
add(const ChainRange& chain)
{
    // TODO: skip the computation entirely if we already know this is positive (in case of the clearing optimization)
    //std::list<Entry> c(std::begin(chain), std::end(chain));
    //c.sort([this](const Entry& e1, const Entry& e2) { return this->cmp_(e1.index(), e2.index()); });
    auto entry_cmp = [this](const Entry& e1, const Entry& e2) { return this->cmp_(e1.index(), e2.index()); };
    std::set<Entry,decltype(entry_cmp)>     c(std::begin(chain), std::end(chain), entry_cmp);

    visitors_chain_initialized(c);

    Index pair = Reduction<Index>::reduce(c, reduced_, pairs_, field_,
                                          [this](FieldElement m, Index cl)
                                          { this->visitors_addto<>(m, cl); },
                                          cmp_);

    if (pair != unpaired)
        pairs_[pair] = pairs_.size();

    pairs_.push_back(pair);
    reduced_.emplace_back(std::begin(c), std::end(c));
    visitors_reduction_finished<>();

    return pair;
}

template<class F, typename I, class C, template<class Self> class... V>
template<class ChainRange>
typename dionysus::ReducedMatrix<F,I,C,V...>::Index
dionysus::ReducedMatrix<F,I,C,V...>::
set(Index i, const ChainRange& chain)
{
    std::list<Entry> c(std::begin(chain), std::end(chain));
    c.sort([this](const Entry& e1, const Entry& e2) { return this->cmp_(e1.index(), e2.index()); });
    visitors_chain_initialized(c);

    Index pair = Reduction<Index>::reduce(c, reduced_, pairs_, field_,
                                          [this](FieldElement m, Index cl)
                                          { this->visitors_addto<>(m, cl); },
                                          cmp_);

    if (pair != unpaired)
        pairs_[pair] = i;

    pairs_[i]   = pair;
    reduced_[i] = Chain(std::begin(c), std::end(c));
    visitors_reduction_finished<>();

    return pair;
}
