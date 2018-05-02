#include <iostream>
#include <vector>
#include <fstream>

#include <dionysus/simplex.h>
#include <dionysus/filtration.h>
#include <dionysus/fields/zp.h>
#include <dionysus/fields/z2.h>
#include <dionysus/ordinary-persistence.h>
#include <dionysus/cohomology-persistence.h>
#include <dionysus/zigzag-persistence.h>
#include <dionysus/standard-reduction.h>
#include <dionysus/row-reduction.h>
namespace d = dionysus;

#include <grid/grid.h>

#include <cnpy.h>

#include <format.h>
#include <opts/opts.h>

#include "grid-topology.h"

template<unsigned D, class K_ = d::Z2Field>
struct ExecuteEP
{
    typedef         K_                                              K;

    typedef         GridTopology<D>                                 Topology;
    typedef         typename Topology::Vertex                       Vertex;
    typedef         grid::GridRef<float,D>                          Function;
    typedef         typename Function::Value                        Value;
    typedef         std::tuple<Value, Vertex>                       ValueVertex;
    typedef         std::vector<ValueVertex>                        Vertices;

    typedef         d::Simplex<Vertex>                              Simplex;
    typedef         d::Filtration<Simplex>                          Filtration;

    //typedef         d::OrdinaryPersistence<K>                       Persistence;
    //typedef         d::OrdinaryPersistenceNoNegative<K>             Persistence;
    //typedef         d::CohomologyPersistence<K>                     Persistence;
    //typedef         d::ZigzagPersistence<K>                         Persistence;

    //typedef         d::StandardReduction<Persistence>               Reduction;
    typedef         d::RowReduction<K>                              Reduction;


    void            operator()(cnpy::NpyArray& arr, K k, std::ostream* out, bool extended = true) const
    {
        // load the data
        Function        function(reinterpret_cast<Value*>(arr.data), typename Function::Vertex(arr.shape), !arr.fortran_order);
        Topology        topology(function.shape());

        // generate the filtration
        Vertices vertices;
        for(const Vertex& v : topology.vertices())
            vertices.push_back(ValueVertex(function(v), v));
        std::sort(vertices.begin(), vertices.end());

        std::vector<size_t> ops;

        Vertex w = std::numeric_limits<Vertex>::max();
        Filtration filtration { Simplex {w} };
        for(auto& vval : vertices)
        {
            Value   val;
            Vertex  v;
            std::tie(val,v) = vval;

            std::vector<Simplex>    lower_star = topology.lower_star(v, function);
            std::sort(lower_star.begin(), lower_star.end());        // order by dimension

            for (auto&& s : lower_star)
                filtration.emplace_back(s);

            ops.push_back(filtration.size());
        }

        if (extended)
            for (auto it = vertices.rbegin(); it != vertices.rend(); ++it)
            {
                Vertex v = std::get<1>(*it);
                std::vector<Simplex>    upper_star = topology.upper_star(v, function);
                std::sort(upper_star.begin(), upper_star.end());

                for (auto& s : upper_star)
                    filtration.emplace_back(s.join(w));

                ops.push_back(filtration.size());
            }

        fmt::print("Filtration size: {}\n", filtration.size());

        // compute persistence
        typedef         typename Reduction::Index               Index;

        //Persistence                     persistence(k);
        //Reduction                       reduce(persistence);
        Reduction                       reduce(k);
        reduce(filtration,
               [out,&ops,&vertices](int dimension, Index birth, Index death)
               {
                   size_t birth_v = std::upper_bound(ops.begin(), ops.end(), birth) - ops.begin();
                   size_t death_v = std::upper_bound(ops.begin(), ops.end(), death) - ops.begin();
                   if (birth_v == death_v)
                       return;

                   Value birth_val;
                   if (birth_v < vertices.size())
                       birth_val = std::get<0>(vertices[birth_v]);
                   else
                       birth_val = std::get<0>(vertices[2*vertices.size() - birth_v]);

                   Value death_val;
                   if (death_v < vertices.size())
                       death_val = std::get<0>(vertices[death_v]);
                   else
                       death_val = std::get<0>(vertices[2*vertices.size() - death_v - 1]);

                   fmt::print(*out, "{} {} {}\n",
                               dimension, birth_val, death_val);
               });
        fmt::print("Reduction finished\n");
    }
};

int main(int argc, char** argv)
{
    using namespace opts;
    Options ops;

    bool no_extended, help;
    ops
        >> Option('n', "no-extended", no_extended, "don't compute the extended part")
        >> Option('h', "help",        help,        "show help")
    ;

    std::string infn;
    if (!ops.parse(argc,argv) || help || !(ops >> PosOption(infn)))
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
    if (arr.word_size != sizeof(ExecuteEP<0>::Value))
    {
        fmt::print("Word sizes don't match: {} vs {}\n", arr.word_size, sizeof(ExecuteEP<0>::Value));
        return 1;
    }

    if (arr.shape.size() == 3)
        ExecuteEP<3>()(arr, d::Z2Field(), out, !no_extended);
    else if (arr.shape.size() == 2)
        ExecuteEP<2>()(arr, d::Z2Field(), out, !no_extended);
    else
        fmt::print("Can't process array of dimension: {}\n", arr.shape.size());

    arr.destruct();
}
