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
        typedef         typename RowMatrix::IndexPair               IndexPair;
        typedef         typename RowMatrix::FieldElement            FieldElement;
        typedef         typename RowMatrix::IndexChain              IndexChain;
        typedef         typename RowMatrix::Column                  Column;

        typedef         std::unordered_map<Index, IndexChain>       Matrix;
        typedef         std::unordered_map<Index, Index>            BirthIndexMap;


                        ZigzagPersistence(const Field&      field,
                                          const Comparison& cmp = Comparison()):
                            Z(field, cmp), B(field, cmp),
                            cell_indices(0),
                            z_indicies_last(0),
                            z_indicies_first(-1),
                            b_indices(0)                            {}

        template<class ChainRange>
        Index           add(const ChainRange& chain);               // returns the id of the dying cycle (or unpaired)
        //Index           remove(Index i);

        void            reserve(size_t)                             {}              // here for compatibility only
        const Field&    field() const                               { return Z.field(); }
        const Comparison&   cmp() const                             { return Z.cmp(); }

        static
        const Index     unpaired = Reduction<Index>::unpaired;

    private:
        RowMatrix       Z, B;
        Matrix          D;

        BirthIndexMap   birth_index;
        Index           cell_indices;
        Index           z_indicies_last, z_indicies_first;
        Index           b_indices;
};

}

#include "zigzag-persistence.hpp"

#endif
