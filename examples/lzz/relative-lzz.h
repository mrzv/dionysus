#ifndef DIONYSUS_RELATIVE_LZZ_H
#define DIONYSUS_RELATIVE_LZZ_H

#include <vector>
#include <tuple>

#include <boost/range/adaptors.hpp>
namespace ba = boost::adaptors;

#include <dionysus/relative-homology-zigzag.h>
#include <dionysus/simplex.h>
namespace d = dionysus;

template<class Field_, class Topology_, class Function_>
struct RelativeLZZ
{
    typedef                 Field_                                              Field;
    typedef                 Topology_                                           Topology;
    typedef                 Function_                                           Function;

    typedef                 d::RelativeHomologyZigzag<Field>                    Zigzag;
    typedef                 typename Zigzag::Index                              Index;

    typedef                 typename Topology::Vertex                           Vertex;
    typedef                 typename Topology::Simplex                          Simplex;
    typedef                 d::ChainEntry<Field, Simplex>                       CellChainEntry;
    typedef                 d::ChainEntry<Field, Index>                         ChainEntry;

    typedef                 typename Function::Value                            Value;
    typedef                 std::tuple<Value, Vertex>                           ValueVertex;
    typedef                 std::vector<ValueVertex>                            Vertices;

    typedef                 std::tuple<int, Index>                              CofacesIndex;
    typedef                 std::unordered_map<Simplex, CofacesIndex>           Complex;

                            RelativeLZZ(const Field& field_, const Topology& topology_, const Function& function_):
                                topology(topology_), function(function_), zz(field_)    {}

    void                    add_both(const Simplex& s);
    void                    remove_both(const Simplex& s);

    Index                   remove_relative(const Simplex& s);
    Index                   add_relative(const Simplex& s);

    template<class ReportPair>
    void                    operator()(const ReportPair& report_pair);

    const Topology&         topology;
    const Function&         function;

    Zigzag                  zz;

    Complex                 cmplx;
    Index                   idx = 0;

    Index                   op = 0;
    std::vector<Index>      ops;            // stores a list of the last Index added for the given vertex (to convert indices to vertices)

    std::unordered_map<Index,bool>      birth_type_;        // records whether index gives birth on addition or removal
};

#include "relative-lzz.hpp"

#endif
