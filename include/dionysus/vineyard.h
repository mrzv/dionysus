#pragma once

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "chain.h"
#include "trails-chains.h"

namespace dionysus
{

template<class Field_, typename Index_ = unsigned>
class Vineyard;

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

        friend class Vineyard<Field_, Index_>;

        Field       field_;
        Chains      columns_;
        Indices     cell_at_;
        Indices     position_;
        Indices     low_of_column_;
        Indices     column_with_low_;
};

template<class Field_, typename Index_>
class Vineyard
{
    public:
        using Field        = Field_;
        using Index        = Index_;
        using FieldElement = typename Field::Element;
        using Entry        = ChainEntry<Field, Index>;
        using Chain        = std::vector<Entry>;
        using Chains       = std::vector<Chain>;
        using Indices      = std::vector<Index>;
        using LocalCells   = std::array<Index, 4>;
        using Matrix       = VineyardMatrix<Field, Index>;
        using Persistence  = OrdinaryPersistenceWithV<Field, Index>;

    public:
        Vineyard(const Field& field, Chains boundary):
            field_(field),
            reduced_(field_, boundary.size()),
            chains_(),
            pairs_(boundary.size(), unpaired())
        {
            for (Index i = 0; i < size(); ++i)
                validate_chain(boundary[i]);

            Persistence persistence(field_);
            persistence.reserve(size());
            for (Index i = 0; i < size(); ++i)
                persistence.add(std::move(boundary[i]));

            chains_ = std::move(persistence.template visitor<0>().v_);

            for (Index i = 0; i < size(); ++i)
            {
                reduced_.columns_[i] = std::move(persistence.column(i));
                reduced_.set_low_unchecked(i, initial_low(reduced_.columns_[i]));
                pairs_[i] = persistence.pair(i);
            }
        }

        size_t          size() const                            { return reduced_.size(); }
        const Field&    field() const                           { return field_; }

        const Chain&    reduced_column(Index column) const       { return reduced_[column]; }
        const Chain&    chain(Index column) const                { return chains_.at(column); }

        Index           cell_at(Index p) const                   { return reduced_.cell_at(p); }
        Index           position(Index cell) const               { return reduced_.position(cell); }
        Index           low(Index column) const                  { return reduced_.low(column); }
        Index           pivot(Index row) const                   { return reduced_.pivot(row); }

        Index           pair(Index i) const                      { return pairs_.at(i); }

        std::pair<Index, Index> transpose_position(Index p)
        {
            if (p + 1 >= size())
                throw std::out_of_range("vineyard transposition position out of range");

            Index a = reduced_.cell_at(p);
            Index b = reduced_.cell_at(p + 1);
            LocalCells local = {{ a, b, pairs_[a], pairs_[b] }};
            bool a_positive = reduced_.low(a) == unpaired();
            bool b_positive = reduced_.low(b) == unpaired();
            Index low_a = reduced_.low(a);
            Index low_b = reduced_.low(b);
            Index pivot_a = reduced_.pivot(a);
            Index pivot_b = reduced_.pivot(b);
            bool r_pivot_b_contains_a = pivot_b != unpaired() && contains(reduced_.columns_[pivot_b], a);

            bool cancelled_v = cancel_chain_entry(b, a);
            auto swapped = reduced_.transpose_position(p);

            if (a_positive && b_positive)
            {
                if (pivot_a != unpaired() && pivot_b != unpaired() && r_pivot_b_contains_a)
                {
                    if (reduced_.position(pivot_a) < reduced_.position(pivot_b))
                        add_to_cancel_low(pivot_b, pivot_a, a);
                    else
                        add_to_cancel_low(pivot_a, pivot_b, a);
                }
            } else if (!a_positive && !b_positive)
            {
                if (cancelled_v && reduced_.position(low_b) < reduced_.position(low_a))
                    add_to_cancel_low(a, b, low_a);
            } else if (!a_positive && b_positive)
            {
                if (cancelled_v)
                    add_to_cancel_low(a, b, low_a);
            }

            refresh_local_lows(local);
            refresh_local_pairs(local);
            return swapped;
        }

        static Index    unpaired()                              { return Matrix::unpaired(); }

    private:
        void validate_index(Index i) const
        {
            if (i >= size())
                throw std::out_of_range("vineyard index out of range");
        }

        void validate_chain(const Chain& chain) const
        {
            for (const Entry& e : chain)
                validate_index(e.index());
        }

        void refresh_local_lows(const LocalCells& local)
        {
            for (size_t i = 0; i < local.size(); ++i)
            {
                Index column = local[i];
                if (skip_local(local, i, column))
                    continue;
                reduced_.set_low_unchecked(column, reduced_.compute_low(reduced_.columns_[column]));
            }
        }

        void refresh_local_pairs(const LocalCells& local)
        {
            for (size_t i = 0; i < local.size(); ++i)
            {
                Index cell = local[i];
                if (skip_local(local, i, cell))
                    continue;

                Index partner = pairs_[cell];
                if (partner != unpaired())
                    pairs_[partner] = unpaired();
                pairs_[cell] = unpaired();
            }

            for (size_t i = 0; i < local.size(); ++i)
            {
                Index column = local[i];
                if (skip_local(local, i, column))
                    continue;

                Index l = reduced_.low(column);
                if (l == unpaired())
                    continue;
                if (!contains(local, l))
                    throw std::logic_error("vineyard local pair update escaped local set");
                if ((pairs_[column] != unpaired() && pairs_[column] != l) ||
                    (pairs_[l] != unpaired() && pairs_[l] != column))
                    throw std::logic_error("vineyard pair update is inconsistent");
                pairs_[column] = l;
                pairs_[l] = column;
            }
        }

        bool skip_local(const LocalCells& local, size_t i, Index cell) const
        {
            if (cell == unpaired())
                return true;
            validate_index(cell);
            for (size_t j = 0; j < i; ++j)
                if (local[j] == cell)
                    return true;
            return false;
        }

        bool contains(const LocalCells& local, Index cell) const
        {
            for (Index local_cell : local)
                if (local_cell == cell)
                    return true;
            return false;
        }

        void add_to_cancel_low(Index column, Index other, Index low)
        {
            FieldElement x = coefficient(reduced_.columns_[column], low);
            FieldElement y = coefficient(reduced_.columns_[other], low);
            FieldElement m = field_.neg(field_.div(x, y));
            add_to(column, m, other);
        }

        bool cancel_chain_entry(Index column, Index other)
        {
            FieldElement x;
            if (!coefficient(chains_[column], other, x))
                return false;

            FieldElement y = coefficient(chains_[other], other);
            FieldElement m = field_.neg(field_.div(x, y));
            add_to(column, m, other);
            return true;
        }

        void add_to(Index column, FieldElement m, Index other)
        {
            dionysus::Chain<Chain>::addto(reduced_.columns_[column], m, reduced_.columns_[other], field_, entry_cmp);
            dionysus::Chain<Chain>::addto(chains_[column], m, chains_[other], field_, entry_cmp);
        }

        FieldElement coefficient(const Chain& chain, Index row) const
        {
            FieldElement element;
            if (coefficient(chain, row, element))
                return element;
            throw std::logic_error("vineyard low entry missing from column");
        }

        bool coefficient(const Chain& chain, Index row, FieldElement& element) const
        {
            auto it = std::lower_bound(chain.begin(), chain.end(), row,
                                       [](const Entry& e, Index i) { return e.index() < i; });
            if (it == chain.end() || it->index() != row)
                return false;
            element = it->element();
            return true;
        }

        bool contains(const Chain& chain, Index row) const
        {
            auto it = std::lower_bound(chain.begin(), chain.end(), row,
                                       [](const Entry& e, Index i) { return e.index() < i; });
            return it != chain.end() && it->index() == row;
        }

        static bool entry_cmp(const Entry& x, const Entry& y)     { return x.index() < y.index(); }

        static Index initial_low(const Chain& chain)
        {
            if (chain.empty())
                return unpaired();
            return chain.back().index();
        }

        Field       field_;
        Matrix      reduced_;
        Chains      chains_;
        Indices     pairs_;
};

}
