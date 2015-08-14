#ifndef DIONYSUS_GRID_TOPOLOGY_H
#define DIONYSUS_GRID_TOPOLOGY_H

#include <boost/range/counting_range.hpp>

#include <dionysus/simplex.h>
#include <dionysus/rips.h>          // for bron-kebosch to generate the star in the Freudenthal triangulation
namespace d = dionysus;

#include <grid/box.h>
#include <grid/vertices.h>

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

#endif
