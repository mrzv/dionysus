#pragma once

#include <cassert>
#include <vector>
#include <iostream>

namespace dionysus
{


template<class Matrix_>
class MatrixFiltrationCell;


// adapt Matrix as a Filtration to make it possible to feed into reduction algorithms
template <class Matrix_>
class MatrixFiltration
{
    public:
        using Matrix = Matrix_;
        using Dimensions = std::vector<short unsigned>;
        using Cell = MatrixFiltrationCell<Matrix>;

    public:
                MatrixFiltration(const Matrix* m, Dimensions dimensions):
                    m_(m), dimensions_(dimensions)              { assert(m_->size() == dimensions_.size()); }

        Cell            operator[](size_t i) const      { return Cell(this, i); }
        size_t          size() const                    { return m_->size(); }

    private:
        const Matrix*   m_;
        Dimensions      dimensions_;

        friend class MatrixFiltrationCell<Matrix>;
};


template<class Matrix_>
class MatrixFiltrationCell
{
    public:
        using Matrix = Matrix_;
        using Field = typename Matrix::Field;
        using MatrixFiltration = MatrixFiltration<Matrix>;
        using ChainEntry = ChainEntry<Field,MatrixFiltrationCell>;
        using BoundaryChain = std::vector<ChainEntry>;

    public:
                MatrixFiltrationCell(const MatrixFiltration* mf, size_t i):
                    mf_(mf), i_(i)      {}

        short unsigned  dimension() const       { return mf_->dimensions_[i_]; }

        bool            operator==(const MatrixFiltrationCell& other) const     { return i_ == other.i_; }
        bool            operator!=(const MatrixFiltrationCell& other) const     { return i_ != other.i_; }

        BoundaryChain   boundary() const
        {
            BoundaryChain bdry;
            for (auto& entry : (*mf_->m_)[i_])
                bdry.emplace_back(ChainEntry { entry.e, MatrixFiltrationCell(mf_, entry.i) });
            return bdry;
        }

        size_t          i() const       { return i_; }

        friend
        std::ostream&   operator<<(std::ostream& out, const MatrixFiltrationCell& c)
        { out << c.i_; return out; }

    private:
        const MatrixFiltration* mf_ = nullptr;
        size_t                  i_;
};


}
