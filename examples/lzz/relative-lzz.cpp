#include <vector>
#include <iostream>
#include <fstream>

#include <boost/range/counting_range.hpp>
#include <boost/algorithm/minmax.hpp>

#include <dionysus/fields/zp.h>
#include <dionysus/fields/z2.h>
#include <dionysus/rips.h>          // for bron-kebosch to generate the star in the Freudenthal triangulation
namespace d = dionysus;

#include "relative-lzz.h"

#include <grid/grid.h>
#include <grid/box.h>
#include <grid/vertices.h>

#include <cnpy.h>

#include <format.h>
#include <opts/opts.h>

#if 0
struct DummyTopology
{
    typedef     int                         Vertex;
    typedef     d::Simplex<>                Simplex;

    std::vector<Vertex>         vertices() const                { return { 0, 1, 2 }; }
    std::vector<Simplex>        closed_star(Vertex v) const     { return { {0}, {1}, {2}, {0,1}, {0,2}, {1,2}, {0,1,2} }; }
};

struct DummyFunction
{
    typedef     float                       Value;

    Value       operator()(int v) const     { return v; }
};
#endif

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

                GridTopology(const Position& shape);

    Vertices    vertices() const                                { return Vertices(Vertex(0), box.position_to_vertex()(box.to()) + 1); }
    Simplices   closed_star(Vertex v) const;

    Box                 box;
    PositionSimplices   star;
};

template<unsigned D>
GridTopology<D>::
GridTopology(const Position& shape): box(shape)
{
    struct FakeDistances
    {
        typedef     Position        IndexType;
        typedef     float           DistanceType;       // meaningless
    };
    typedef     d::Rips<FakeDistances,PositionSimplex>          Rips;
    typedef     typename Rips::VertexContainer                  VertexContainer;

    VertexContainer vertices(grid::VerticesIterator<Position>::begin(-Position::one(), Position::one()),
                             grid::VerticesIterator<Position>::end  (-Position::one(), Position::one()));

    auto neighbor = [](const Position& x, const Position& y)        // neighbor
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
                    };

    VertexContainer current;
    Rips::bron_kerbosch(current, vertices, std::prev(vertices.begin()), D,
                        neighbor,
                        [this,neighbor](PositionSimplex&& s)
                        {
                            const PositionSimplex& cs = s;
                            for (auto& u : cs)
                                if (!neighbor(u,Position::zero()))
                                    return;
                            star.emplace_back(s);
                        });

    //fmt::print("Star generated; size = {}\n", star.size());
    //for (auto& s : star)
    //    fmt::print("   {}\n", s);
}

template<unsigned D>
typename GridTopology<D>::Simplices
GridTopology<D>::
closed_star(Vertex v) const
{
    Simplices   result;
    Position    p = box.position(v);
    auto        bounds = box.bounds_test();
    auto        vertex = box.position_to_vertex();
    for (auto& ps : star)
    {
        bool accept = true;
        for (auto& u : ps)
            if (!bounds(p + u))
            {
                accept = false;
                break;
            }

        if (!accept) continue;

        result.emplace_back(std::make_pair(ps.begin(),ps.end()) | ba::transformed([vertex,p](const Position& u) { return vertex(p+u); }));
    }
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
