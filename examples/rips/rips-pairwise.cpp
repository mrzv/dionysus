#include <iostream>
#include <vector>
#include <fstream>

#include <dionysus/distances.h>
#include <dionysus/rips.h>
#include <dionysus/filtration.h>
#include <dionysus/fields/zp.h>
#include <dionysus/fields/z2.h>
#include <dionysus/ordinary-persistence.h>
#include <dionysus/cohomology-persistence.h>
#include <dionysus/zigzag-persistence.h>
#include <dionysus/standard-reduction.h>
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

//typedef         d::Z2Field                                              K;
typedef         d::ZpField<>                                            K;
//typedef         d::OrdinaryPersistence<K>                               Persistence;
typedef         d::CohomologyPersistence<K>                             Persistence;
//typedef         d::ZigzagPersistence<K>                                 Persistence;

int main(int argc, char* argv[])
{
    using opts::Options;
    using opts::Option;
    using opts::Present;
    using opts::PosOption;

    short unsigned          skeleton = 2;
    DistanceType            max_distance = std::numeric_limits<DistanceType>::infinity();
    std::string             infilename, diagram_name;

    Options ops(argc, argv);
    ops
        >> Option('s', "skeleton",      skeleton,           "dimension of the Rips complex we want to compute")
        >> Option('m', "max-distance",  max_distance,       "maximum distance value cutoff")
    ;

    if ( (ops >> Present('h', "help", "show help message") ||
        !(ops >> PosOption(infilename))))
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
    std::cout << "Generated complex of size: " << filtration.size() << std::endl;

    // Generate filtration with respect to distance and compute its persistence
    filtration.sort(Generator::Comparison(distances));

    //K k;
    K k(11);
    Persistence                             persistence(k);
    d::StandardReduction<Persistence>       reduce(persistence);
    reduce(filtration);
    std::cout << "Reduction finished" << std::endl;
}
