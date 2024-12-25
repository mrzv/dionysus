#ifndef HERA_AUCTION_PARAMS_H
#define HERA_AUCTION_PARAMS_H

#include <unordered_map>
#include <limits>
#include <cassert>
#include <iostream>
#include <ostream>
#include <ios>

#include <hera/common.h>

namespace hera {

template<class Real_ = double>
struct AuctionParams {
    using Real = Real_;
    Real wasserstein_power {1};
    Real delta {0.01}; // relative error
    Real internal_p {get_infinity<Real>()};
    Real initial_epsilon {0}; // 0.0 means maxVal / 4.0
    Real epsilon_common_ratio {5};
    int max_num_phases {std::numeric_limits<decltype(max_num_phases)>::max()};
    int max_bids_per_round {1};  // imitate Gauss-Seidel is default behaviour
    unsigned int dim {2}; // for pure geometric version only; ignored in persistence diagrams
    bool tolerate_max_iter_exceeded {false}; // whether auction should throw an exception on max. iterations exceeded
    bool return_matching {false}; // whether to return optimal matching along with cost
    bool match_inf_points {true}; // whether to add infinite points to matching; ignored, if return_matching is false
};

template<class Real>
std::ostream& operator<<(std::ostream& out, const AuctionParams<Real>& p)
{
    out << "AuctionParams(dim=" << p.dim << ", wasserstein_power=" << p.wasserstein_power << ", delta=" << p.delta << ", internal_p=";
    if (is_infinity(p.internal_p))
        out << "INF";
    else
        out << p.internal_p;
    out << ", initial_epsilon=" << p.initial_epsilon << ", epsilon_common_ratio=" << p.epsilon_common_ratio;
    out << ", max_num_phases=" << p.max_num_phases << ", max_bids_per_round=" << p.max_bids_per_round;
    out << std::boolalpha;
    out << ", tolerate_max_iter_exceeded=" << p.tolerate_max_iter_exceeded;
    out << ", return_matching=" << p.return_matching;
    out << ", match_inf_points=" << p.match_inf_points << ")";

    return out;
}


} // namespace hera

#endif //HERA_AUCTION_PARAMS_H
