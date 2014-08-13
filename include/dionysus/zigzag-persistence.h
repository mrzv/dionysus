#ifndef DIONYSUS_ZIGZAG_PERSISTENCE_H
#define DIONYSUS_ZIGZAG_PERSISTENCE_H

#include <tuple>

#include "sparse-row-matrix.h"

namespace dionysus
{

template<class Field_, class Index_ = int, class Comparison_ = std::less<Index_>>
class ZigzagPersistence
{
    public:
        typedef         Field_                                      Field;
        typedef         Index_                                      Index;
        typedef         Comparison_                                 Comparison;

        typedef         SparseRowMatrix<Field, Index, Comparison>   RowMatrix;
        typedef         SparseRowMatrix<Field, Index, Comparison,
                                        std::deque>                 DequeRowMatrix;
        typedef         typename RowMatrix::IndexPair               IndexPair;
        typedef         typename RowMatrix::FieldElement            FieldElement;
        typedef         typename RowMatrix::IndexChain              IndexChain;
        typedef         typename RowMatrix::Column                  Column;
        typedef         typename RowMatrix::Row                     Row;
        typedef         typename DequeRowMatrix::Column             DequeColumn;

        typedef         std::unordered_map<Index, Index>            BirthIndexMap;


                        ZigzagPersistence(const Field&      field,
                                          const Comparison& cmp = Comparison()):
                            Z(field, cmp), B(field, cmp), D(field, cmp),
                            operations(0),
                            cell_indices(0),
                            z_indicies_last(0),
                            z_indicies_first(-1),
                            b_indices(0)                            {}

        template<class ChainRange>
        Index           add(const ChainRange& chain);               // returns the id of the dying cycle (or unpaired)
        Index           remove(Index cell);

        void            reserve(size_t)                             {}              // here for compatibility only
        const Field&    field() const                               { return Z.field(); }
        const Comparison&   cmp() const                             { return Z.cmp(); }

        static
        const Index     unpaired = Reduction<Index>::unpaired;

    private:
        RowMatrix       Z, D;
        DequeRowMatrix  B;

        BirthIndexMap   birth_index;
        Index           operations;
        Index           cell_indices;
        Index           z_indicies_last, z_indicies_first;
        Index           b_indices;
};

}

#include "zigzag-persistence.hpp"

#endif
