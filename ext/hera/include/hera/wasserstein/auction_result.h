#ifndef HERA_AUCTION_RESULT_H
#define HERA_AUCTION_RESULT_H

#include <unordered_map>
#include <cmath>
#include <ostream>

#include <hera/common.h>

namespace hera {

template<class Real = double>
struct AuctionResult {
    Real cost {0};                  // cost of the final matching
    Real distance {0};              // Wasserstein distance, cost^(1/q)
    Real final_relative_error {0};  // the real relative error, can be less than user-specified delta
    Real start_epsilon {-1};        // epsilon used in the first phase of epsilon-scaling
    Real final_epsilon {-1};        // epsilon used in the last phase of epsilon-scaling
    long int num_rounds {0};        // total number of bidding rounds in all phases
    int num_phases {0};             // number of epsilon-scaling phases
    std::vector<Real> prices;       // final prices

    void compute_distance(Real q)  { distance = std::pow(cost, 1/ q); }

    std::unordered_map<int, int> matching_a_to_b_;
    std::unordered_map<int, int> matching_b_to_a_;

    void clear_matching()
    {
        matching_a_to_b_.clear();
        matching_b_to_a_.clear();
    }

    void add_to_matching(int a, int b)
    {
        assert(matching_a_to_b_.count(a) == 0 and matching_b_to_a_.count(b) == 0);
        matching_a_to_b_[a] = b;
        matching_b_to_a_[b] = a;
    }
};


template<class Real>
std::ostream& operator<<(std::ostream& out, const AuctionResult<Real>& r)
{
    out << "Result(cost=" << r.cost << ", distance=" << r.distance << ", num_rounds=" << r.num_rounds << ", num_phases=" << r.num_phases << ", error=" << r.final_relative_error;
    out << ", start_epsilon=" << r.start_epsilon << ", final_epsilon=" << r.final_epsilon << ")";
    return out;
}


template<class Real>
AuctionResult<Real> add_results(const AuctionResult<Real>& r1, const AuctionResult<Real>& r2, Real q)
{
    AuctionResult<Real> result;
    result.cost = r1.cost + r2.cost;
    result.compute_distance(q);
    result.num_rounds = r1.num_rounds + r2.num_rounds;
    result.num_phases = r1.num_phases + r2.num_phases;

    result.matching_a_to_b_ = r1.matching_a_to_b_;
    result.matching_a_to_b_.insert(r2.matching_a_to_b_.begin(), r2.matching_a_to_b_.end());

    result.matching_b_to_a_ = r1.matching_b_to_a_;
    result.matching_b_to_a_.insert(r2.matching_b_to_a_.begin(), r2.matching_b_to_a_.end());

    // we add only results from matching infinite points, where relative error is 0, to finite
    // TODO: fix this, now it is an upper bound
    result.final_relative_error = std::max(result.final_relative_error, r2.final_relative_error);

    return result;
}

} // namespace hera

#endif //HERA_AUCTION_RESULT_H
