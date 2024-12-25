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
#include <locale>
#include <iomanip>
#include <vector>

#include <opts/opts.h>

#include <hera/wasserstein.h>

// any container of pairs of Reals can be used,
// we use vector in this example.

int main(int argc, char* argv[])
{
    using Real = double;
    using PairVector = std::vector<std::pair<Real, Real>>;
    PairVector diagramA, diagramB;

    hera::AuctionParams<Real> params;
    params.max_num_phases = 800;
    bool help { false };
    bool print_relative_error { false };

    opts::Options ops;
    ops >> opts::Option('q', "degree", params.wasserstein_power, "Wasserstein degree")
        >> opts::Option('d', "error", params.delta, "Relative error")
        >> opts::Option('p', "internal-p", params.internal_p, "Internal norm")
        >> opts::Option("initial-epsilon", params.initial_epsilon, "Initial epsilon")
        >> opts::Option("epsilon-factor", params.epsilon_common_ratio, "Epsilon factor")
        >> opts::Option("max-bids-per-round", params.max_bids_per_round, "Maximal number of bids per round")
        >> opts::Option('m', "max-rounds", params.max_num_phases, "Maximal number of iterations")
        >> opts::Option('e', "--print-error", print_relative_error, "Print real relative error")
        >> opts::Option('t', "tolerate", params.tolerate_max_iter_exceeded, "Suppress max-iterations-exceeded error and print the best result.")
        >> opts::Option('h', "help", help, "Print help")
    ;


    std::string dgm_fname_1, dgm_fname_2;
    if (!ops.parse(argc, argv) || help || !(ops >> opts::PosOption(dgm_fname_1) >> opts::PosOption(dgm_fname_2))) {
        std::cerr << "Usage: " << argv[0] << " file-1 file-2\n" << std::endl;
        std::cerr << "compute Wasserstein distance between persistence diagrams in file1 and file2.\n";
        std::cout << ops << std::endl;
        return 1;
    }

    if (!hera::read_diagram_point_set<Real, PairVector>(dgm_fname_1, diagramA)) {
        std::exit(1);
    }

    if (!hera::read_diagram_point_set(dgm_fname_2, diagramB)) {
        std::exit(1);
    }

    if (params.wasserstein_power < 1.0) {
        std::cerr << "Wasserstein_degree was \"" << params.wasserstein_power << "\", must be a number >= 1.0. Cannot proceed. " << std::endl;
        std::exit(1);
    }

    if (params.wasserstein_power == 1.0) {
        hera::remove_duplicates<Real>(diagramA, diagramB);
    }

    if ( params.delta <= 0.0) {
        std::cerr << "relative error was \"" << params.delta << "\", must be a number > 0.0. Cannot proceed. " << std::endl;
        std::exit(1);
    }

    // default for internal metric is l_infinity
    if (std::isinf(params.internal_p)) {
        params.internal_p = hera::get_infinity<Real>();
    }

    if (not hera::is_p_valid_norm<Real>(params.internal_p)) {
        std::cerr << "internal-p was \"" << params.internal_p << "\", must be a number >= 1.0 or inf. Cannot proceed. " << std::endl;
        std::exit(1);
    }

    // if you want to specify initial value for epsilon and the factor
    // for epsilon-scaling
    if (params.initial_epsilon < 0.0) {
        std::cerr << "initial-epsilon was \"" << params.initial_epsilon << "\", must be a non-negative number. Cannot proceed." << std::endl;
        std::exit(1);
    }

    if (params.epsilon_common_ratio <= 1.0 and params.epsilon_common_ratio != 0.0) {
        std::cerr << "The 7th argument (epsilon factor) was \"" << params.epsilon_common_ratio << "\", must be a number greater than 1. Cannot proceed." << std::endl;
        std::exit(1);
    }

    if (params.max_bids_per_round == 0)
        params.max_bids_per_round = std::numeric_limits<decltype(params.max_bids_per_round)>::max();

    auto res = hera::wasserstein_cost_detailed(diagramA, diagramB, params);

    std::cout << std::setprecision(15) << res.distance << std::endl;
    if (print_relative_error)
        std::cout << "Relative error: " << res.final_relative_error << std::endl;

    return 0;

}
