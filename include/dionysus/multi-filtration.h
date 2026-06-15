#ifndef DIONYSUS_MULTI_FILTRATION_H
#define DIONYSUS_MULTI_FILTRATION_H

#include <vector>

#include "multi-filtration-common.h"

namespace b   = boost;
namespace bmi = boost::multi_index;

namespace dionysus
{

// MultiFiltration stores a filtered cell complex with duplicate-preserving lookup.
// It allows for bidirectional translation between a cell and its ordered index.
template<class Cell_,
         bool  checked_index = false>
class MultiFiltration:
    private detail::DuplicateFiltrationStorage<detail::DuplicateCellWithIndex<Cell_>, checked_index>
{
    private:
        using Storage = detail::DuplicateFiltrationStorage<detail::DuplicateCellWithIndex<Cell_>, checked_index>;

    public:
        using               order = typename Storage::order;

        typedef             Cell_                                               Cell;
        using               CellWithIndex = detail::DuplicateCellWithIndex<Cell>;
        using               CellLookupIndex = typename Storage::CellLookupIndex;
        using               Container = typename Storage::Container;
        typedef             typename Storage::value_type                        value_type;

        typedef             typename Storage::Complex                           Complex;
        typedef             typename Storage::Order                             Order;

        using               OrderCell = typename Storage::OrderCell;
        using               OrderConstIterator = typename Storage::OrderConstIterator;
        using               OrderIterator = typename Storage::OrderIterator;

    public:
                            MultiFiltration()                                        = default;
                            MultiFiltration(MultiFiltration&& other)                 = default;
        MultiFiltration&    operator=(MultiFiltration&& other)                       = default;

                            MultiFiltration(const std::initializer_list<Cell>& cells):
                                MultiFiltration(std::begin(cells), std::end(cells))  {}

        template<class Iterator>
                            MultiFiltration(Iterator bg, Iterator end)               { for (auto it = bg; it != end; ++it) push_back(*it); }

        template<class CellRange>
                            MultiFiltration(const CellRange& cells):
                                MultiFiltration(std::begin(cells), std::end(cells))  {}

        // Lookup
        const Cell&         operator[](size_t i) const                                { return Storage::operator[](i); }
        size_t              index(const Cell& s, size_t i) const                     { return this->preceding_index(s, i); }
        bool                contains(const Cell& s) const                            { return Storage::contains(s); }

        void                push_back(const Cell& s)                                 { this->push_indexed(CellWithIndex(s, this->next_index())); }
        void                push_back(Cell&& s)                                      { this->push_indexed(CellWithIndex(std::move(s), this->next_index())); }

        void                replace(size_t i, const Cell& s)                         { this->replace_indexed(i, CellWithIndex(s, i)); }

        template<class... Args>
        void                emplace_back(Args&&... args)                             { this->push_indexed(CellWithIndex(Cell(std::forward<Args>(args)...), this->next_index())); }

        template<class Cmp = std::less<Cell>>
        void                sort(const Cmp& cmp = Cmp())                             { this->sort_with_updates(cmp, [](CellWithIndex& c, size_t i) { c.i = i; }); }

        void                rearrange(const std::vector<size_t>& indices)            { this->rearrange_with_updates(indices, [](CellWithIndex& c, size_t i) { c.i = i; }); }

        OrderConstIterator  begin() const                                            { return Storage::begin(); }
        OrderConstIterator  end() const                                              { return Storage::end(); }
        OrderIterator       begin()                                                  { return Storage::begin(); }
        OrderIterator       end()                                                    { return Storage::end(); }
        size_t              size() const                                             { return Storage::size(); }
        void                clear()                                                  { Storage::clear(); }

        Cell&               back()                                                   { return Storage::back(); }
        const Cell&         back() const                                             { return Storage::back(); }
};

}

#endif
