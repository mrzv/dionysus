/*
 
Copyright (c) 2015, M. Kerber, D. Morozov, A. Nigmetov
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
(Enhancements) to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to copyright holder,
without imposing a separate written license agreement for such Enhancements,
then you hereby grant the following license: a  non-exclusive, royalty-free
perpetual license to install, use, modify, prepare derivative works, incorporate
into other computer software, distribute, and sublicense such enhancements or
derivative works thereof, in binary and source code form.

  */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <random>

#include "wasserstein.h"

// any container of pairs of doubles can be used,
// we use vector in this example.

int main(int argc, char* argv[])
{
    geom_ws::PairVector diagramA, diagramB;

    if (argc < 3 ) {
        std::cerr << "Usage: " << argv[0] << " file1 file2 [wasserstein_degree] [relative_error] [internal norm] [output_actual_error]. By default power is 1.0, relative error is 0.01, internal norm is l_infinity, actual relative error is not printed." << std::endl;
        return 1;
    }

    if (!geom_ws::readDiagramPointSet(argv[1], diagramA)) {
        std::exit(1);
    }

    if (!geom_ws::readDiagramPointSet(argv[2], diagramB)) {
        std::exit(1);
    }

    double wasserPower = (4 <= argc) ? atof(argv[3]) : 1.0;
    if (wasserPower < 1.0) {
        std::cerr << "The third argument (wasserstein_degree) was \"" << argv[3] << "\", must be a number >= 1.0. Cannot proceed. " << std::endl;
        std::exit(1);
    }

    if (wasserPower == 1.0) {
        geom_ws::removeDuplicates(diagramA, diagramB);
    }

    //default relative error:  1%
    double delta = (5 <= argc) ? atof(argv[4]) : 0.01;
    if ( delta <= 0.0) {
        std::cerr << "The 4th argument (relative error) was \"" << argv[4] << "\", must be a number > 0.0. Cannot proceed. " << std::endl;
        std::exit(1);
    }

    // default for internal metric is l_infinity
    double internal_p = ( 6 <= argc ) ? atof(argv[5]) : std::numeric_limits<double>::infinity();
    if (internal_p < 1.0) {
        std::cerr << "The 5th argument (internal norm) was \"" << argv[5] << "\", must be a number >= 1.0. Cannot proceed. " << std::endl;
        std::exit(1);
    }

    // if you want to specify initial value for epsilon and the factor
    // for epsilon-scaling
    double initialEpsilon= ( 7 <= argc ) ? atof(argv[6]) : 0.0 ;
    double epsFactor = ( 8 <= argc ) ? atof(argv[7]) : 0.0 ;

    double res = geom_ws::wassersteinDist(diagramA, diagramB, wasserPower, delta, internal_p, initialEpsilon, epsFactor);
    std::cout << std::setprecision(15) << res << std::endl;
    return 0;
}
