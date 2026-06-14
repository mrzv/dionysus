#ifndef DIONYSUS_BOUNDARY_MATRIX_H
#define DIONYSUS_BOUNDARY_MATRIX_H

#include <utility>

namespace dionysus
{

namespace detail
{

template<class MatrixEntry>
typename MatrixEntry::Element signed_boundary_element(typename MatrixEntry::Element element, short prime)
{
    if (element > prime / 2)
        element -= prime;
    return element;
}

}

template<class MatrixFiltration, class Filtration>
MatrixFiltration make_boundary_matrix_filtration(const Filtration& f, short prime = 3)
{
    using Matrix = typename MatrixFiltration::Matrix;
    using Entry = typename Matrix::Entry;
    using Index = typename Matrix::Index;

    Matrix m(prime);
    typename MatrixFiltration::Dimensions dimensions;
    typename MatrixFiltration::Values values;

    m.resize(f.size());
    dimensions.resize(f.size());
    values.resize(f.size());

    Index i = 0;
    for (auto& c : f)
    {
        dimensions[i] = c.dimension();
        values[i] = c.data();

        typename Matrix::Chain column;
        for (auto e : c.boundary(m.field()))
        {
            auto element = detail::signed_boundary_element<Entry>(e.element(), prime);
            column.emplace_back(element, f.index(e.index(), i));
        }
        m.set(i, std::move(column));

        ++i;
    }

    return MatrixFiltration(std::move(m), std::move(dimensions), std::move(values));
}

template<class MatrixFiltration, class Filtration>
MatrixFiltration make_coboundary_matrix_filtration(const Filtration& f, short prime = 3)
{
    using Matrix = typename MatrixFiltration::Matrix;
    using Entry = typename Matrix::Entry;
    using Index = typename Matrix::Index;

    Matrix m(prime);
    typename MatrixFiltration::Dimensions dimensions;
    typename MatrixFiltration::Values values;

    Index n = static_cast<Index>(f.size());
    m.resize(n);
    dimensions.resize(n);
    values.resize(n);

    Index i = 0;
    for (auto& c : f)
    {
        dimensions[n - 1 - i] = c.dimension();
        values[n - 1 - i] = c.data();

        for (auto e : c.boundary(m.field()))
        {
            auto element = detail::signed_boundary_element<Entry>(e.element(), prime);
            auto index = static_cast<Index>(f.index(e.index(), i));
            m.column(n - 1 - index).emplace_back(element, n - 1 - i);
        }

        ++i;
    }

    for (Index i = 0; i < m.size(); ++i)
        m.sort(m.column(i));

    return MatrixFiltration(std::move(m), std::move(dimensions), std::move(values));
}

}

#endif
