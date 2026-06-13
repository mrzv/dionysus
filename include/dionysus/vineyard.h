#pragma once

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "chain.h"
#include "trails-chains.h"

namespace dionysus
{

template<class Field_, typename Index_, class Decomposition_>
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

        std::pair<Index, Index> transpose(Index p)
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

        template<class, typename, class>
        friend class Vineyard;

        Field       field_;
        Chains      columns_;
        Indices     cell_at_;
        Indices     position_;
        Indices     low_of_column_;
        Indices     column_with_low_;
};

struct VineyardVDecomposition {};
struct VineyardUDecomposition {};

template<class Decomposition, class Field, typename Index>
struct VineyardPersistence;

template<class Field, typename Index>
struct VineyardPersistence<VineyardVDecomposition, Field, Index>
{
    using Persistence = OrdinaryPersistenceWithV<Field, Index>;

    static typename Persistence::Chains take_basis(Persistence& persistence)
    {
        return std::move(persistence.template visitor<0>().v_);
    }
};

template<class Field, typename Index>
struct VineyardPersistence<VineyardUDecomposition, Field, Index>
{
    using Persistence = OrdinaryPersistenceWithU<Field, Index>;

    static typename Persistence::Chains take_basis(Persistence& persistence)
    {
        return std::move(persistence.template visitor<0>().u_);
    }
};

template<class Field_, typename Index_, class Decomposition_>
class Vineyard
{
    public:
        using Field        = Field_;
        using Index        = Index_;
        using Decomposition = Decomposition_;
        using FieldElement = typename Field::Element;
        using Entry        = ChainEntry<Field, Index>;
        using Chain        = std::vector<Entry>;
        using Chains       = std::vector<Chain>;
        using Indices      = std::vector<Index>;
        using Matrix       = VineyardMatrix<Field, Index>;
        using PersistenceTraits = VineyardPersistence<Decomposition, Field, Index>;
        using Persistence  = typename PersistenceTraits::Persistence;

    public:
        Vineyard(const Field& field, Chains boundary, bool lazy = false):
            field_(field),
            reduced_(field_, boundary.size()),
            basis_(),
            pairs_(boundary.size(), unpaired()),
            lazy_(lazy)
        {
            for (Index i = 0; i < size(); ++i)
                validate_chain(boundary[i]);

            Persistence persistence(field_);
            persistence.reserve(size());
            for (Index i = 0; i < size(); ++i)
                persistence.add(std::move(boundary[i]));

            basis_ = PersistenceTraits::take_basis(persistence);

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
        const Chain&    basis(Index column) const                { return basis_.at(column); }

        Index           cell_at(Index p) const                   { return reduced_.cell_at(p); }
        Index           position(Index cell) const               { return reduced_.position(cell); }
        Index           low(Index column) const                  { return reduced_.low(column); }
        Index           pivot(Index row) const                   { return reduced_.pivot(row); }

        Index           pair(Index i) const                      { return pairs_.at(i); }

        std::pair<Index, Index> transpose(Index p)
        {
            if (p + 1 >= size())
                throw std::out_of_range("vineyard transposition position out of range");

            Index a = reduced_.cell_at(p);
            Index b = reduced_.cell_at(p + 1);
            bool a_positive = reduced_.low(a) == unpaired();
            bool b_positive = reduced_.low(b) == unpaired();
            Index low_a = reduced_.low(a);
            Index low_b = reduced_.low(b);
            Index pivot_a = reduced_.pivot(a);
            Index pivot_b = reduced_.pivot(b);
            bool r_pivot_b_contains_a = pivot_b != unpaired() && contains(reduced_.columns_[pivot_b], a);

            bool cancelled = clear_transposition_entry(b, a, Decomposition());
            auto swapped = reduced_.transpose(p);

            if (a_positive && b_positive)
            {
                bool pairing_switched = false;
                if (r_pivot_b_contains_a)
                {
                    if (pivot_a == unpaired())
                    {
                        reduced_.set_low_unchecked(pivot_b, a);
                        pair_cells(a, pivot_b);
                        clear_pair(b);
                        pairing_switched = true;
                    } else if (reduced_.position(pivot_a) < reduced_.position(pivot_b))
                    {
                        add_to_cancel_low(pivot_b, pivot_a, a);
                        reduced_.set_low_unchecked(pivot_a, a);
                        reduced_.set_low_unchecked(pivot_b, b);
                    }
                    else
                    {
                        add_to_cancel_low(pivot_a, pivot_b, a);
                        reduced_.set_low_unchecked(pivot_b, a);
                        reduced_.set_low_unchecked(pivot_a, b);
                        pair_cells(a, pivot_b);
                        pair_cells(b, pivot_a);
                        pairing_switched = true;
                    }
                }

                if (!pairing_switched && pivot_a != unpaired() && pivot_b != unpaired() &&
                    reduced_.position(pivot_b) < reduced_.position(pivot_a))
                {
                    if (lazy_)
                        clear_lazy_entry(pivot_a, pivot_b, Decomposition());
                }
            } else if (!a_positive && !b_positive)
            {
                if (cancelled && reduced_.position(low_b) < reduced_.position(low_a))
                {
                    add_to_cancel_low(a, b, low_a);
                    reduced_.set_low_unchecked(a, low_b);
                    reduced_.set_low_unchecked(b, low_a);
                    pair_cells(a, low_b);
                    pair_cells(b, low_a);
                } else
                {
                    reduced_.set_low_unchecked(a, low_a);
                    reduced_.set_low_unchecked(b, low_b);
                    if (lazy_)
                        clear_lazy_entry(b, a, Decomposition());
                }
            } else if (!a_positive && b_positive)
            {
                if (cancelled)
                {
                    add_to_cancel_low(a, b, low_a);
                    reduced_.set_low_unchecked(a, unpaired());
                    reduced_.set_low_unchecked(b, low_a);
                    if (pivot_b != unpaired())
                        reduced_.set_low_unchecked(pivot_b, a);
                    pair_cells(b, low_a);
                    if (pivot_b != unpaired())
                        pair_cells(a, pivot_b);
                    else
                        clear_pair(a);
                }
            }

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

        void clear_pair(Index cell)
        {
            if (cell == unpaired())
                return;
            validate_index(cell);

            Index partner = pairs_[cell];
            if (partner != unpaired())
            {
                validate_index(partner);
                if (pairs_[partner] != cell)
                    throw std::logic_error("vineyard pair map is inconsistent");
                pairs_[partner] = unpaired();
            }
            pairs_[cell] = unpaired();
        }

        void pair_cells(Index x, Index y)
        {
            if (x == unpaired() || y == unpaired())
                throw std::logic_error("vineyard cannot pair an unpaired cell");
            if (x == y)
                throw std::logic_error("vineyard cannot pair a cell with itself");
            validate_index(x);
            validate_index(y);

            clear_pair(x);
            clear_pair(y);
            pairs_[x] = y;
            pairs_[y] = x;
        }

        void add_to_cancel_low(Index column, Index other, Index low)
        {
            FieldElement x = coefficient(reduced_.columns_[column], low);
            FieldElement y = coefficient(reduced_.columns_[other], low);
            FieldElement m = field_.neg(field_.div(x, y));
            add_to(column, m, other);
        }

        bool clear_transposition_entry(Index column, Index other, VineyardVDecomposition)
        {
            FieldElement x;
            if (!coefficient(basis_[column], other, x))
                return false;

            FieldElement y = coefficient(basis_[other], other);
            FieldElement m = field_.neg(field_.div(x, y));
            add_to(column, m, other);
            return true;
        }

        bool clear_transposition_entry(Index column, Index other, VineyardUDecomposition)
        {
            FieldElement x;
            if (!coefficient(basis_[other], column, x))
                return false;

            FieldElement y = coefficient(basis_[column], column);
            FieldElement m = field_.neg(field_.div(x, y));
            add_to(column, field_.neg(m), other);
            return true;
        }

        bool clear_lazy_entry(Index column, Index other, VineyardVDecomposition)
        {
            FieldElement x;
            if (!coefficient(basis_[column], other, x))
                return false;

            FieldElement y = coefficient(basis_[other], other);
            FieldElement m = field_.neg(field_.div(x, y));
            add_to(column, m, other);
            return true;
        }

        bool clear_lazy_entry(Index column, Index other, VineyardUDecomposition)
        {
            FieldElement x;
            if (!coefficient(basis_[other], column, x))
                return false;

            FieldElement y = coefficient(basis_[column], column);
            FieldElement m = field_.div(x, y);
            add_to(column, m, other);
            return true;
        }

        void add_to(Index column, FieldElement m, Index other)
        {
            dionysus::Chain<Chain>::addto(reduced_.columns_[column], m, reduced_.columns_[other], field_, entry_cmp);
            add_to_basis(column, m, other, Decomposition());
        }

        void add_to_basis(Index column, FieldElement m, Index other, VineyardVDecomposition)
        {
            dionysus::Chain<Chain>::addto(basis_[column], m, basis_[other], field_, entry_cmp);
        }

        void add_to_basis(Index column, FieldElement m, Index other, VineyardUDecomposition)
        {
            dionysus::Chain<Chain>::addto(basis_[other], field_.neg(m), basis_[column], field_, entry_cmp);
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
        Chains      basis_;
        Indices     pairs_;
        bool        lazy_;
};

template<class Field_, typename Index_ = unsigned>
using VineyardV = Vineyard<Field_, Index_, VineyardVDecomposition>;

template<class Field_, typename Index_ = unsigned>
using VineyardU = Vineyard<Field_, Index_, VineyardUDecomposition>;

}
