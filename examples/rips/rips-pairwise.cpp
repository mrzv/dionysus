#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <limits>

#include <dionysus/distances.h>
#include <dionysus/rips.h>
#include <dionysus/filtration.h>
#include <dionysus/fields/zp.h>
#include <dionysus/fields/z2.h>
#include <dionysus/ordinary-persistence.h>
#include <dionysus/cohomology-persistence.h>
#include <dionysus/zigzag-persistence.h>
#include <dionysus/standard-reduction.h>
#include <dionysus/row-reduction.h>
#include <dionysus/pair-recorder.h>
namespace d = dionysus;

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
typedef         d::Filtration<Simplex>                                  Filtration;

typedef         d::Z2Field                                              K;
//typedef         d::ZpField<>                                            K;
//typedef         d::OrdinaryPersistence<K>                               Persistence;
typedef         d::PairRecorder<d::CohomologyPersistence<K>>            Persistence;
//typedef         d::ZigzagPersistence<K>                                 Persistence;

int main(int argc, char* argv[])
{
    using opts::Options;
    using opts::Option;
    using opts::PosOption;

    short unsigned          skeleton = 2;
    DistanceType            max_distance = std::numeric_limits<DistanceType>::infinity();
    std::string             infilename, diagram_name;
    int                     p = -1;
    bool                    output_diagram, verbose, help;

    Options ops;
    ops
        >> Option('s', "skeleton",      skeleton,           "dimension of the Rips complex we want to compute")
        >> Option('m', "max-distance",  max_distance,       "maximum distance value cutoff")
        >> Option('p', "p",             p,                  "restrict diagrams to this dimension")
        >> Option('d', "diagram",       output_diagram,     "output diagram")
        >> Option('v', "verbose",       verbose,            "verbose output")
    ;

    if (!ops.parse(argc,argv) || help || !(ops >> PosOption(infilename)))
    {
        std::cout << "Usage: " << argv[0] << " input-points" << std::endl;
        std::cout << ops;
        return 1;
    }

    PointContainer          points;
    read_points(infilename, points);

    PairDistances           distances(points);

    Generator               rips(distances);
    Generator::Evaluator    size(distances);
    Filtration              filtration;

    rips.generate(skeleton, max_distance,
                  [&filtration](Simplex&& s) { filtration.push_back(s); });
    if (verbose) std::cout << "# Generated complex of size: " << filtration.size() << std::endl;

    // Generate filtration with respect to distance and compute its persistence
    filtration.sort(Generator::Comparison(distances));

    K k;
    //K k(11);
    Persistence                             persistence(k);
    d::StandardReduction<Persistence>       reduce(persistence);
    //d::RowReduction<K>                      reduce(k);
    //const auto&                             persistence = reduce.persistence();
    reduce(filtration);
    if (verbose) std::cout << "# Reduction finished" << std::endl;

    std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1);

    if (output_diagram)
    {
        typedef decltype(persistence.pair(0)) Index;
        Generator::Evaluator    eval(distances);
        for (Index i = 0; i < persistence.size(); ++i)
        {
            Index j = persistence.pair(i);
            if (j < i) continue;

            if (filtration[i].dimension() == skeleton)
                continue;

            if (p != -1 && filtration[i].dimension() != p)
                continue;

            auto birth = eval(filtration[i]);
            auto death = birth;

            if (j != persistence.unpaired())
                death = eval(filtration[j]);

            if (birth == death)
                continue;

            std::cout << birth << " ";
            if (j == persistence.unpaired())
                std::cout << "inf\n";
            else
                std::cout << death << '\n';
        }
    }
}
