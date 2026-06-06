#pragma once

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "chain.h"

namespace dionysus
{

template<class Field_, typename Index_ = unsigned>
class VineyardMatrix
{
    public:
        using Field        = Field_;
        using Index        = Index_;
        using FieldElement = typename Field::Element;
        using Entry        = ChainEntry<Field, Index>;
        using Chain        = std::vector<Entry>;
        using Chains       = std::vector<Chain>;
        using Indices      = std::vector<Index>;

    public:
        VineyardMatrix(const Field& field, size_t size):
            field_(field),
            columns_(size),
            cell_at_(size),
            position_(size),
            low_of_column_(size, unpaired()),
            column_with_low_(size, unpaired())
        {
            for (Index i = 0; i < size; ++i)
                cell_at_[i] = position_[i] = i;
        }

        size_t          size() const                            { return columns_.size(); }
        const Field&    field() const                           { return field_; }

        const Chain&    operator[](Index column) const           { return columns_.at(column); }
        const Chain&    column(Index column) const               { return columns_.at(column); }

        Index           cell_at(Index p) const                   { return cell_at_.at(p); }
        Index           position(Index cell) const               { return position_.at(cell); }
        Index           low(Index column) const                  { return low_of_column_.at(column); }
        Index           pivot(Index row) const                   { return column_with_low_.at(row); }

        static Index    unpaired()                              { return std::numeric_limits<Index>::max(); }

        void set_column(Index column, Chain chain)
        {
            check_index(column);
            for (const Entry& e : chain)
                check_index(e.index());

            std::sort(chain.begin(), chain.end(), entry_cmp);
            columns_[column] = std::move(chain);
            set_low_unchecked(column, compute_low(columns_[column]));
        }

        bool contains(Index column, Index row) const
        {
            check_index(column);
            if (row >= size())
                return false;

            const Chain& c = columns_[column];
            auto it = std::lower_bound(c.begin(), c.end(), row,
                                       [](const Entry& e, Index i) { return e.index() < i; });
            return it != c.end() && it->index() == row;
        }

        std::pair<Index, Index> transpose_position(Index p)
        {
            if (p + 1 >= size())
                throw std::out_of_range("vineyard transposition position out of range");

            Index a = cell_at_[p];
            Index b = cell_at_[p + 1];
            std::swap(cell_at_[p], cell_at_[p + 1]);
            std::swap(position_[a], position_[b]);
            return std::make_pair(a, b);
        }

        void refresh_low(Index column)
        {
            check_index(column);
            set_low_unchecked(column, compute_low(columns_[column]));
        }

        void set_low(Index column, Index low)
        {
            check_index(column);
            if (low != unpaired())
                check_index(low);

            set_low_unchecked(column, low);
        }

    private:
        void set_low_unchecked(Index column, Index low)
        {
            Index old_low = low_of_column_[column];
            if (old_low != unpaired() && column_with_low_[old_low] == column)
                column_with_low_[old_low] = unpaired();

            low_of_column_[column] = low;
            if (low != unpaired())
                column_with_low_[low] = column;
        }

        void check_index(Index i) const
        {
            if (i >= size())
                throw std::out_of_range("vineyard matrix index out of range");
        }

        Index compute_low(const Chain& c) const
        {
            if (c.empty())
                return unpaired();

            auto low = std::max_element(c.begin(), c.end(), [this](const Entry& x, const Entry& y)
                                        { return position_[x.index()] < position_[y.index()]; });
            return low->index();
        }

        static bool entry_cmp(const Entry& x, const Entry& y)     { return x.index() < y.index(); }
        Field       field_;
        Chains      columns_;
        Indices     cell_at_;
        Indices     position_;
        Indices     low_of_column_;
        Indices     column_with_low_;
};

}
