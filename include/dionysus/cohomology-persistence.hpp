template<class F, class I, class Cmp>
template<class ChainRange>
typename dionysus::CohomologyPersistence<F,I,Cmp>::Index
dionysus::CohomologyPersistence<F,I,Cmp>::
add(const ChainRange& chain)
{
    auto entry_cmp = [this](const Entry& e1, const Entry& e2) { return this->cmp_(e1.index(), e2.index()); };
    std::set<Entry,decltype(entry_cmp)>     row_sum(entry_cmp);
    for (auto it = std::begin(chain); it != std::end(chain); ++it)
        for (auto& re : rows_[it->index()])
            Chain<decltype(row_sum)>::addto(row_sum, it->element(), Entry(re.element(), re.column->index(), re.column), field_, cmp_);

    if (row_sum.empty())        // Birth
    {
        columns_.emplace_back(rows_.size());
        auto before_end = columns_.end();
           --before_end;
        columns_.back().chain.push_back(Entry(field_.id(), rows_.size(), before_end));
        rows_.emplace_back();
        rows_.back().push_back(columns_.back().chain.front());
        return unpaired;
    } else                      // Death
    {
        // Select front element in terms of comparison (rows are unsorted)
        auto it = std::max_element(std::begin(row_sum), std::end(row_sum), entry_cmp);

        Entry first = std::move(*it);
        row_sum.erase(it);

        for (auto& ce : row_sum)
        {
            FieldElement ay = field_.neg(field_.div(ce.element(), first.element()));
            std::list<Entry> tmp(ce.column->chain.begin(), ce.column->chain.end());

            Chain<decltype(tmp)>::addto(tmp, ay, first.column->chain, field_, cmp_);

            ce.column->chain = Column(std::make_move_iterator(std::begin(tmp)),
                                      std::make_move_iterator(std::end(tmp)));
            for (auto it = ce.column->chain.begin(); it != ce.column->chain.end(); ++it)
            {
                it->column = ce.column;
                rows_[it->index()].push_back(*it);
            }
        }
        Index pair = first.column->index();
        columns_.erase(first.column);
        rows_.emplace_back();       // useless row; only present to make indices match
        return pair;
    }
}
