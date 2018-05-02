#include <vector>
#include <iostream>
#include <fstream>

#include <dionysus/fields/zp.h>
#include <dionysus/fields/z2.h>

#include "relative-lzz.h"

#include <grid/grid.h>

#include <cnpy.h>

#include <format.h>
#include <opts/opts.h>

#include "grid-topology.h"


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
    Options ops;

    bool help;
    std::string infn;

    ops >> Option('h', "help", help, "show help");

    if (ops.parse(argc,argv) || help || !(ops >> PosOption(infn)))
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
