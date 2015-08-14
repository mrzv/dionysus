#include <vector>
#include <iostream>
#include <fstream>

#include <boost/range/counting_range.hpp>
#include <boost/algorithm/minmax.hpp>

#include <dionysus/fields/zp.h>
#include <dionysus/fields/z2.h>
#include <dionysus/rips.h>          // for bron-kebosch to generate the star in the Freudenthal triangulation
namespace d = dionysus;

#include <format.h>
#include "relative-lzz.h"

#include <grid/grid.h>
#include <grid/box.h>
#include <grid/vertices.h>

#include <cnpy.h>

#include <format.h>
#include <opts/opts.h>

template<unsigned D>
struct GridTopology
{
    typedef     grid::Box<D>                                    Box;
    typedef     typename Box::Vertex                            Vertex;
    typedef     typename Box::Position                          Position;
    typedef     decltype(boost::counting_range(Vertex(0), Vertex(1)))   Vertices;
    typedef     d::Simplex<Position>                            PositionSimplex;
    typedef     std::vector<PositionSimplex>                    PositionSimplices;
    typedef     d::Simplex<Vertex>                              Simplex;
    typedef     std::vector<Simplex>                            Simplices;

    struct FakeDistances
    {
        typedef     Position        IndexType;
        typedef     float           DistanceType;       // meaningless
    };
    typedef     d::Rips<FakeDistances,Simplex>                  Rips;
    typedef     typename Rips::VertexContainer                  VertexContainer;


                GridTopology(const Position& shape);

    Vertices    vertices() const                                { return Vertices(Vertex(0), box.position_to_vertex()(box.to()) + 1); }

    template<class Function>
    Simplices   neighborhood(Vertex v, const Function& f, bool lower, bool include) const;

    template<class Function>
    Simplices   upper_star(Vertex v, const Function& f) const   { return neighborhood(v,f,false,true); }
    template<class Function>
    Simplices   upper_link(Vertex v, const Function& f) const   { return neighborhood(v,f,false,false); }
    template<class Function>
    Simplices   lower_star(Vertex v, const Function& f) const   { return neighborhood(v,f,true,true); }
    template<class Function>
    Simplices   lower_link(Vertex v, const Function& f) const   { return neighborhood(v,f,true,false); }

    static bool neighbor(const Position& x, const Position& y);
    bool        neighbor_vertex(Vertex x, Vertex y) const       { return neighbor(box.position(x), box.position(y)); }

    Box                     box;
    std::vector<Position>   neighbors;
};

template<unsigned D>
GridTopology<D>::
GridTopology(const Position& shape): box(shape)
{
    auto vi  = grid::VerticesIterator<Position>::begin(-Position::one(), Position::one()),
         end = grid::VerticesIterator<Position>::end  (-Position::one(), Position::one());
    while (vi != end)
    {
        if (neighbor(*vi,Position::zero()) && *vi != Position::zero())
            neighbors.push_back(*vi);
        ++vi;
    }
}

template<unsigned D>
bool
GridTopology<D>::
neighbor(const Position& x, const Position& y)
{
    Position diff = x - y;

    int min = diff[0], max = diff[0];
    for (unsigned i = 0; i < D; ++i)
    {
        if (diff[i] < min)
            min = diff[i];
        if (diff[i] > max)
            max = diff[i];
    }

    if (min < -1) return false;
    if (max > +1) return false;

    return !(min < 0 && max > 0);
}

template<unsigned D>
template<class Function>
typename GridTopology<D>::Simplices
GridTopology<D>::
neighborhood(Vertex v, const Function& f, bool lower, bool include) const
{
    typedef                 typename Function::Value                            Value;
    typedef                 std::tuple<Value, Vertex>                           ValueVertex;

    Position    vp = box.position(v);
    ValueVertex vval(f(v), v);

    auto vertex = box.position_to_vertex();

    VertexContainer vertices;
    for (auto up : neighbors)
    {
        up += vp;
        if (!box.contains(up))
            continue;
        Vertex u = vertex(up);
        ValueVertex uval(f(up), u);
        if ((lower && uval < vval) || (!lower && uval > vval))
            vertices.push_back(u);
    }

    Simplices       result;
    VertexContainer current;
    if (include)
        current.push_back(v);
    Rips::bron_kerbosch(current, vertices, std::prev(vertices.begin()), D,
                        [this](Vertex u, Vertex v) { return neighbor_vertex(u,v); },
                        [this,&result](Simplex&& s) { result.emplace_back(s); });

    return result;
}


template<unsigned D, class K_ = d::Z2Field>
struct ExecuteLZZ
{
    typedef         K_                                              K;

    typedef         GridTopology<D>                                 Topology;
    typedef         grid::GridRef<float,D>                          Function;
    typedef         typename Function::Value                        Value;
    typedef         typename Function::Vertex                       Vertex;

    typedef         RelativeLZZ<K, Topology, Function>              RLZZ;
    typedef         typename RLZZ::ValueVertex                      ValueVertex;

    void            operator()(cnpy::NpyArray& arr, K k, std::ostream* out) const
    {
        // load the data
        Function        function(reinterpret_cast<Value*>(arr.data), Vertex(arr.shape), !arr.fortran_order);
        Topology        topology(function.shape());

        RLZZ            lzz(k, topology, function);
        lzz([out](int dimension, const ValueVertex& birth, const ValueVertex& death, bool birth_type, bool death_type)
            {
                fmt::print(*out, "{}{} {} {} {}\n",
                            birth_type ? '+' : '-',
                            death_type ? '+' : '-',
                            dimension, std::get<0>(birth), std::get<0>(death));
            });
    }
};

int main(int argc, char** argv)
{
    using namespace opts;
    Options ops(argc, argv);

    std::string infn;
    if (  ops >> Present('h', "help", "show help") ||
        !(ops >> PosOption(infn)))
    {
        fmt::print("Usage: {} IN.npy [OUT.dgm]\n{}", argv[0], ops);
        return 1;
    }

    std::string outfn;
    std::ostream* out = &std::cout;
    std::ofstream ofs;
    if (ops >> PosOption(outfn))
    {
        ofs.open(outfn.c_str());
        out = &ofs;
    }

    cnpy::NpyArray arr = cnpy::npy_load(infn);
    if (arr.word_size != sizeof(ExecuteLZZ<0>::Value))
    {
        fmt::print("Word sizes don't match: {} vs {}\n", arr.word_size, sizeof(ExecuteLZZ<0>::Value));
        return 1;
    }

    if (arr.shape.size() == 3)
        ExecuteLZZ<3>()(arr, d::Z2Field(), out);
    else if (arr.shape.size() == 2)
        ExecuteLZZ<2>()(arr, d::Z2Field(), out);
    else
        fmt::print("Can't process array of dimension: {}\n", arr.shape.size());

    arr.destruct();
}
