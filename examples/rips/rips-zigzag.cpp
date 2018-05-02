#include <iostream>
#include <vector>
#include <unordered_map>

#include <boost/range/adaptors.hpp>
namespace ba = boost::adaptors;

#include <dionysus/simplex.h>
#include <dionysus/fields/zp.h>
#include <dionysus/fields/z2.h>
#include <dionysus/distances.h>
#include <dionysus/rips.h>
#include <dionysus/zigzag-persistence.h>
namespace d = dionysus;

#include <dionysus/dlog/progress.h>

#include <opts/opts.h>

#include <common.h>     // read_points()

typedef         std::vector<float>                                      Point;
typedef         std::vector<Point>                                      PointContainer;

typedef         d::PairwiseDistances<PointContainer,
                                     d::L2Distance<Point>>              PairDistances;
typedef         PairDistances::DistanceType                             DistanceType;
typedef         PairDistances::IndexType                                Vertex;

typedef         d::Rips<PairDistances>                                  Generator;
typedef         Generator::Simplex                                      Simplex;
typedef         std::set<Simplex>                                       SimplexSet;

typedef         std::vector<Vertex>                                     VertexVector;
typedef         std::vector<DistanceType>                               EpsilonVector;
typedef         std::tuple<Vertex,Vertex>                               Edge;
typedef         std::vector<Edge>                                       EdgeVector;

//typedef         d::Z2Field                                              K;
typedef         d::ZpField<>                                            K;
typedef         d::Simplex<>                                            Simplex;
typedef         d::ZigzagPersistence<K>                                 Persistence;
typedef         typename Persistence::Index                             Index;

typedef         std::unordered_map<Simplex, Index>                      Complex;
typedef         d::ChainEntry<K, Simplex>                               SimplexChainEntry;
typedef         d::ChainEntry<K, Index>                                 ChainEntry;

// debug
typedef         std::unordered_map<Index, Simplex>                      RComplex;

// Information we need to know when a class dies
struct      BirthInfo
{
    typedef         short unsigned                                      Dimension;

                    BirthInfo(DistanceType dist = DistanceType(), Dimension dim = Dimension()):
                        distance(dist), dimension(dim)              {}
    DistanceType    distance;
    Dimension       dimension;
};

typedef         std::unordered_map<Index, BirthInfo>                    BirthMap;


int main(int argc, char** argv)
{
    using opts::Options;
    using opts::Option;
    using opts::PosOption;

    short unsigned          skeleton = 2;
    DistanceType            multiplier = 6;
    short unsigned          p = 11;
    std::string             infilename, diagram_name;
    bool                    help;

    Options ops;
    ops
        >> Option('s', "skeleton",      skeleton,           "dimension of the Rips complex we want to compute")
        >> Option('m', "multiplier",    multiplier,         "multiplier for epsilon (distance to the next maxmin point)")
        >> Option('p', "prime",         p,                  "prime for arithmetic")
        >> Option('h', "help",          help,               "show help message")
    ;

    if (!ops.parse(argc,argv) || !(ops >> PosOption(infilename)) || !(ops >> PosOption(diagram_name)))
    {
        std::cout << "Usage: " << argv[0] << " input-points diagram.out" << std::endl;
        std::cout << ops;
        return 1;
    }

    PointContainer          points;
    read_points(infilename, points);

    std::ofstream   dgm_out(diagram_name);
    std::ostream&   out = dgm_out;

    // Construct distances and Rips generator
    PairDistances           distances(points);
    Generator               rips(distances);
    Generator::Evaluator    size(distances);

    // Order vertices and epsilons (in maxmin fashion)
    VertexVector        vertices;
    EpsilonVector       epsilons;
    EdgeVector          edges;
    DistanceType        inf     = std::numeric_limits<DistanceType>::infinity();

    {
        EpsilonVector   dist(distances.size(), inf);

        vertices.push_back(distances.begin());
        //epsilons.push_back(inf);
        while (vertices.size() < distances.size())
        {
            for (Vertex v = distances.begin(); v != distances.end(); ++v)
                dist[v] = std::min(dist[v], distances(v, vertices.back()));
            auto max = std::max_element(dist.begin(), dist.end());
            vertices.push_back(max - dist.begin());
            epsilons.push_back(*max);
        }
        epsilons.push_back(0);
    }

    // Generate and sort all the edges
    for (unsigned i = 0; i != vertices.size(); ++i)
        for (unsigned j = i+1; j != vertices.size(); ++j)
        {
            Vertex u = vertices[i];
            Vertex v = vertices[j];
            if (distances(u,v) <= multiplier*epsilons[j-1])
                edges.emplace_back(u,v);
        }
    std::sort(edges.begin(), edges.end(),
              [&distances](const Edge& e1, const Edge& e2)
              { return distances(std::get<0>(e1), std::get<1>(e1)) < distances(std::get<0>(e2), std::get<1>(e2)); });

    // Construct zigzag
    //K               k;
    K               k(p);
    Persistence     persistence(k);
    Complex         simplices;
#ifdef DIONYSUS_ZIGZAG_DEBUG
    RComplex        rsimplices;
#endif

    // Insert vertices
    Index       op   = 0;
    Index       cell = 0;
    BirthMap    births;
    for (auto v : vertices)
    {
        // Add a vertex
        Simplex s = {v};

        // We don't actually need to transform the boundary here,
        // since it's empty anyway, but we keep it for the sake of completeness
        Index pair = persistence.add(s.boundary(persistence.field()) |
                                                ba::transformed([&simplices](const SimplexChainEntry& e)
                                                { return ChainEntry(e.element(), simplices.find(e.index())->second); }));

#ifdef DIONYSUS_ZIGZAG_DEBUG
        rsimplices.emplace(cell, s);
        persistence.check_boundaries([&simplices](const Simplex& s) { return simplices[s]; },
                                     [&rsimplices](Index i)         { return rsimplices.find(i)->second; });
#endif

        births[op++] = BirthInfo(0,0);                  // record the birth
        simplices.emplace(std::move(s), cell++);        // record the cell id
    }

    // Process vertices
    dlog::progress progress(vertices.size());
    unsigned    ce = 0;         // index of the current one past last edge in the complex
    SimplexSet  cofaces;        // record the cofaces of all the simplices that need to be removed and reinserted
    for (unsigned stage = 0; stage != vertices.size() - 1; ++stage)
    {
        unsigned i = vertices.size() - 1 - stage;

        /* Increase epsilon */
        cofaces.clear();

        // Add anything else that needs to be inserted into the complex
        while (ce < edges.size())
        {
            Vertex u,v;
            std::tie(u,v) = edges[ce];
            if (distances(u,v) <= multiplier*epsilons[i-1])
                ++ce;
            else
                break;
            //std::cout << "Adding cofaces of " << u << ' ' << v << std::endl;
            rips.edge_cofaces(u, v,
                              skeleton,
                              multiplier*epsilons[i-1],
                              [&cofaces](Simplex&& s) { cofaces.insert(s); },
                              vertices.begin(),
                              vertices.begin() + i + 1);
        }

        // Insert all the cofaces
        for (auto& s : cofaces)
        {
            //std::cout << "Inserting: " << s << std::endl;

            Index pair = persistence.add(s.boundary(persistence.field()) |
                                                    ba::transformed([&simplices](const SimplexChainEntry& e)
                                                    { return ChainEntry(e.element(), simplices.find(e.index())->second); }));
            simplices.emplace(std::move(s), cell);      // record the cell id
#ifdef DIONYSUS_ZIGZAG_DEBUG
            rsimplices.emplace(cell, s);
            persistence.check_boundaries([&simplices](const Simplex& s) { return simplices[s]; },
                                         [&rsimplices](Index i)         { return rsimplices.find(i)->second; });
#endif
            ++cell;

            if (pair == Persistence::unpaired())
                births[op++] = BirthInfo(epsilons[i-1],s.dimension());              // record the birth
            else
            {
                const BirthInfo& birth = births[pair];
                if ((birth.distance - epsilons[i-1]) != 0 && birth.dimension < skeleton)
                    out << birth.dimension << " " << birth.distance << " " << epsilons[i-1] << std::endl;
                births.erase(pair);
                ++op;
            }
        }

        /* Remove the vertex */
        //std::cout << "Removing vertex: " << vertices[i] << std::endl;
        cofaces.clear();
        rips.vertex_cofaces(vertices[i],
                            skeleton,
                            multiplier*epsilons[i-1],
                            [&cofaces](Simplex&& s) { cofaces.insert(s); },
                            vertices.begin(),
                            vertices.begin() + i + 1);
        //std::cout << "Total cofaces: " << cofaces.size() << std::endl;

        for (auto& s : cofaces | ba::reversed)
        {
            //std::cout << "Removing: " << s << std::endl;
            Complex::const_iterator  it    = simplices.find(s);
            Index                    c     = it->second;
            simplices.erase(it);

            Index pair  = persistence.remove(c);
#ifdef DIONYSUS_ZIGZAG_DEBUG
            rsimplices.erase(c);
            persistence.check_boundaries([&simplices](const Simplex& s) { return simplices[s]; },
                                         [&rsimplices](Index i)         { return rsimplices.find(i)->second; });
#endif

            if (pair == Persistence::unpaired())
                births[op++] = BirthInfo(epsilons[i-1],s.dimension() - 1);          // record the birth
            else
            {
                const BirthInfo& birth = births[pair];
                if ((birth.distance - epsilons[i-1]) != 0 && birth.dimension < skeleton)
                    out << birth.dimension << " " << birth.distance << " " << epsilons[i-1] << std::endl;
                births.erase(pair);
                ++op;
            }
        }

        ++progress;
    }

    // Remove the last vertex
    Index pair = persistence.remove(0);
    simplices.erase((Complex::const_iterator) simplices.begin());     // TODO: add an assertion that the complex has only 1 simplex
#ifdef DIONYSUS_ZIGZAG_DEBUG
    rsimplices.erase(0);
    persistence.check_boundaries([&simplices](const Simplex& s) { return simplices[s]; },
                                 [&rsimplices](Index i)         { return rsimplices.find(i)->second; });
#endif

    const BirthInfo& birth = births[pair];
    out << birth.dimension << " " << birth.distance << " " << epsilons[0] << std::endl;
    ++progress;

    std::cout << "Finished" << std::endl;
}
