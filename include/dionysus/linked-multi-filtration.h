#ifndef DIONYSUS_LINKED_MULTI_FILTRATION_H
#define DIONYSUS_LINKED_MULTI_FILTRATION_H

#include <vector>

#include "multi-filtration-common.h"

namespace b   = boost;
namespace bmi = boost::multi_index;

namespace dionysus
{

// LinkedMultiFiltration stores for each cell a linked index. The linked index is checked first by the index(c,i) search.
template<class Cell_,
         bool  checked_index = false>
class LinkedMultiFiltration:
    private detail::DuplicateFiltrationStorage<detail::LinkedDuplicateCellWithIndex<Cell_>, checked_index>
{
    private:
        using Storage = detail::DuplicateFiltrationStorage<detail::LinkedDuplicateCellWithIndex<Cell_>, checked_index>;

    public:
        using               order = typename Storage::order;

        typedef             Cell_                                               Cell;
        using               LinkedCellWithIndex = detail::LinkedDuplicateCellWithIndex<Cell>;
        using               CellLookupIndex = typename Storage::CellLookupIndex;
        using               Container = typename Storage::Container;
        typedef             typename Storage::value_type                        value_type;

        typedef             typename Storage::Complex                           Complex;
        typedef             typename Storage::Order                             Order;

        using               OrderCell = typename Storage::OrderCell;
        using               OrderConstIterator = typename Storage::OrderConstIterator;
        using               OrderIterator = typename Storage::OrderIterator;

    public:
                            LinkedMultiFiltration()                                        = default;
                            LinkedMultiFiltration(LinkedMultiFiltration&& other)           = default;
        LinkedMultiFiltration&
                            operator=(LinkedMultiFiltration&& other)                       = default;

                            LinkedMultiFiltration(const std::initializer_list<Cell>& cells):
                                LinkedMultiFiltration(std::begin(cells), std::end(cells))  {}

        template<class Iterator>
                            LinkedMultiFiltration(Iterator bg, Iterator end)               { for (auto it = bg; it != end; ++it) push_back(*it, this->next_index()); }

        template<class CellRange>
                            LinkedMultiFiltration(const CellRange& cells):
                                LinkedMultiFiltration(std::begin(cells), std::end(cells))  {}

        // Lookup
        const Cell&         operator[](size_t i) const                                      { return Storage::operator[](i); }
        size_t              index(const Cell& s, size_t i) const;
        bool                contains(const Cell& s) const                                  { return Storage::contains(s); }

        void                push_back(const Cell& s, size_t l)                             { this->push_indexed(LinkedCellWithIndex(s, this->next_index(), l)); }
        void                push_back(Cell&& s, size_t l)                                  { this->push_indexed(LinkedCellWithIndex(std::move(s), this->next_index(), l)); }

        void                replace(size_t i, const Cell& s, size_t l)                     { this->replace_indexed(i, LinkedCellWithIndex(s, i, l)); }

        template<class... Args>
        void                emplace_back(size_t l, Args&&... args)                         { this->push_indexed(LinkedCellWithIndex(Cell(std::forward<Args>(args)...), this->next_index(), l)); }

        template<class Cmp = std::less<Cell>>
        void                sort(const Cmp& cmp = Cmp())                                   { this->sort_with_index_map(cmp, [](LinkedCellWithIndex& c, size_t i, const std::vector<size_t>& indices) { c.i = i; c.linked = indices[c.linked]; }); }

        void                rearrange(const std::vector<size_t>& indices)                  { this->rearrange_with_index_map(indices, [](LinkedCellWithIndex& c, size_t i, const std::vector<size_t>& reverse_indices) { c.i = i; c.linked = reverse_indices[c.linked]; }); }

        OrderConstIterator  begin() const                                                  { return Storage::begin(); }
        OrderConstIterator  end() const                                                    { return Storage::end(); }
        OrderIterator       begin()                                                        { return Storage::begin(); }
        OrderIterator       end()                                                          { return Storage::end(); }
        size_t              size() const                                                   { return Storage::size(); }
        void                clear()                                                        { Storage::clear(); }

        Cell&               back()                                                         { return Storage::back(); }
        const Cell&         back() const                                                   { return Storage::back(); }
};

}

template<class C, bool checked_index>
size_t
dionysus::LinkedMultiFiltration<C,checked_index>::
index(const Cell& s, size_t i) const
{
    if (i < size())
    {
        size_t l = this->get_order()[i].linked;
        if (l < size() && (*this)[l] == s)
            return l;
    }

    return this->preceding_index(s, i);
}

#endif
