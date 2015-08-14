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
        return unpaired();
    } else                      // Death
    {
        // Select front element in terms of comparison (rows are unsorted)
        auto it = std::max_element(std::begin(row_sum), std::end(row_sum), entry_cmp);

        Entry first = std::move(*it);
        row_sum.erase(it);

        for (auto& ce : row_sum)
        {
            FieldElement ay = field_.neg(field_.div(ce.element(), first.element()));
            Chain<Column>::addto(ce.column->chain, ay, first.column->chain, field_,
                                 [this](const Entry& e1, const Entry& e2)
                                 { return this->cmp_(e1.index(), e2.index()); });

            for (auto& x : ce.column->chain)
            {
                x.column = ce.column;
                rows_[x.index()].push_back(x);
            }
        }
        Index pair = first.column->index();
        columns_.erase(first.column);
        rows_.emplace_back();       // useless row; only present to make indices match
        return pair;
    }
}
