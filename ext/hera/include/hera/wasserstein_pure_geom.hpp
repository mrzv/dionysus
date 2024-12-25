#ifndef WASSERSTEIN_PURE_GEOM_HPP
#define WASSERSTEIN_PURE_GEOM_HPP

#define WASSERSTEIN_PURE_GEOM


#include "common/diagram_reader.h"
#include "wasserstein/auction_oracle_kdtree_pure_geom.h"
#include "wasserstein/auction_runner_gs.h"
#include "wasserstein/auction_runner_jac.h"

namespace hera
{
namespace ws
{

    template <class Real>
    using DynamicTraits = typename hera::ws::dnn::DynamicPointTraits<Real>;

    template <class Real>
    using DynamicPoint = typename hera::ws::dnn::DynamicPointTraits<Real>::PointType;

    template <class Real>
    using DynamicPointVector = typename hera::ws::dnn::DynamicPointVector<Real>;

    template <class Real>
    using AuctionRunnerGSR = typename hera::ws::AuctionRunnerGS<Real, hera::ws::AuctionOracleKDTreePureGeom<Real>, hera::ws::dnn::DynamicPointVector<Real>>;

    template <class Real>
    using AuctionRunnerJacR = typename hera::ws::AuctionRunnerJac<Real, hera::ws::AuctionOracleKDTreePureGeom<Real>, hera::ws::dnn::DynamicPointVector<Real>>;



inline AuctionResult<double> wasserstein_cost_detailed(const DynamicPointVector<double>& set_A, const DynamicPointVector<double>& set_B, const AuctionParams<double>& params, const std::vector<double>& prices=std::vector<double>())
{
    using Real = double;

    if (params.wasserstein_power < 1.0) {
        throw std::runtime_error("Bad q in Wasserstein " + std::to_string(params.wasserstein_power));
    }

    if (params.delta < 0.0) {
        throw std::runtime_error("Bad delta in Wasserstein " + std::to_string(params.delta));
    }

    if (params.initial_epsilon < 0.0) {
        throw std::runtime_error("Bad initial epsilon in Wasserstein" + std::to_string(params.initial_epsilon));
    }

    if (params.epsilon_common_ratio < 0.0) {
        throw std::runtime_error("Bad epsilon factor in Wasserstein " + std::to_string(params.epsilon_common_ratio));
    }

    if (set_A.size() != set_B.size()) {
        throw std::runtime_error("Different cardinalities of point clouds: " + std::to_string(set_A.size()) + " != " +  std::to_string(set_B.size()));
    }

    DynamicTraits<Real> traits(params.dim);

    if (params.dim == 1) {
        AuctionResult<Real> result;

        std::vector<std::pair<Real, size_t>> set_A_copy, set_B_copy;
        set_A_copy.reserve(set_A.size());
        set_B_copy.reserve(set_B.size());

        for(size_t i = 0; i < set_A.size(); ++i) {
            set_A_copy.emplace_back(set_A[i][0], i);
        }

        for(size_t i = 0; i < set_B.size(); ++i) {
            set_B_copy.emplace_back(set_B[i][0], i);
        }

        std::sort(set_A_copy.begin(), set_A_copy.end());
        std::sort(set_B_copy.begin(), set_B_copy.end());

        for(size_t i = 0; i < set_A_copy.size(); ++i) {
            auto a = std::get<0>(set_A_copy[i]);
            auto b = std::get<0>(set_B_copy[i]);

            if (params.return_matching) {
                int id_a = std::get<1>(set_A_copy[i]);
                int id_b = std::get<1>(set_B_copy[i]);
                result.add_to_matching(id_a, id_b);
            }

            result.cost += std::pow(std::fabs(a - b), params.wasserstein_power);
        }
        result.distance = std::pow(result.cost, Real(1) / params.wasserstein_power);
        return result;
    } else {

        DynamicPointVector<Real> set_A_copy(set_A);
        DynamicPointVector<Real> set_B_copy(set_B);

        // set point id to the index in vector
        for(size_t i = 0; i < set_A_copy.size(); ++i) {
            traits.id(set_A_copy[i]) = i;
            traits.id(set_B_copy[i]) = i;
        }

        if (params.max_bids_per_round == 1) {
            hera::ws::AuctionRunnerGSR<Real> auction(set_A_copy, set_B_copy, params, prices);
            auction.run_auction();
            return auction.get_result();
        } else {
            hera::ws::AuctionRunnerJacR<Real> auction(set_A_copy, set_B_copy, params, prices);
            auction.run_auction();
            return auction.get_result();
        }
    }



}


inline double wasserstein_cost(const DynamicPointVector<double>& set_A, const DynamicPointVector<double>& set_B, const AuctionParams<double>& params, const std::vector<double>& prices=std::vector<double>())
{
    return wasserstein_cost_detailed(set_A, set_B, params, prices).cost;
}


inline double wasserstein_dist(const DynamicPointVector<double>& set_A, const DynamicPointVector<double>& set_B, const AuctionParams<double>& params)
{
    return std::pow(wasserstein_cost(set_A, set_B, params), 1.0 / params.wasserstein_power);
}

} // ws
} // hera


#endif
