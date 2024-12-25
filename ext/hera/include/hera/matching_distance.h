#pragma once

#include <vector>
#include <limits>
#include <utility>
#include <ostream>
#include <chrono>
#include <tuple>
#include <algorithm>


#include "matching/common_defs.h"
#include "matching/cell_with_value.h"
#include "matching/box.h"
#include "matching/dual_point.h"
#include "matching/dual_box.h"
#include "matching/persistence_module.h"
#include "matching/bifiltration.h"
#include "bottleneck.h"

namespace md {

#ifdef MD_PRINT_HEAT_MAP
    template<class Real>
    using HeatMap = std::map<DualPoint<Real>, Real>;

    template<class Real>
    using HeatMaps = std::map<int, HeatMap<Real>>;
#endif

    enum class BoundStrategy {
        bruteforce,
        local_dual_bound,
        local_dual_bound_refined,
        local_dual_bound_for_each_point,
        local_combined
    };

    enum class TraverseStrategy {
        depth_first,
        breadth_first,
        breadth_first_value,
        upper_bound
    };

    inline std::ostream& operator<<(std::ostream& os, const BoundStrategy& s)
    {
        switch(s) {
            case BoundStrategy::bruteforce :
                os << "bruteforce";
                break;
            case BoundStrategy::local_dual_bound :
                os << "local_grob";
                break;
            case BoundStrategy::local_combined :
                os << "local_combined";
                break;
            case BoundStrategy::local_dual_bound_refined :
                os << "local_refined";
                break;
            case BoundStrategy::local_dual_bound_for_each_point :
                os << "local_for_each_point";
                break;
            default:
                os << "FORGOTTEN BOUND STRATEGY";
        }
        return os;
    }

    inline std::ostream& operator<<(std::ostream& os, const TraverseStrategy& s)
    {
        switch(s) {
            case TraverseStrategy::depth_first :
                os << "DFS";
                break;
            case TraverseStrategy::breadth_first :
                os << "BFS";
                break;
            case TraverseStrategy::breadth_first_value :
                os << "BFS-VAL";
                break;
            case TraverseStrategy::upper_bound :
                os << "UB";
                break;
            default:
                os << "FORGOTTEN TRAVERSE STRATEGY";
        }
        return os;
    }

    inline std::istream& operator>>(std::istream& is, TraverseStrategy& s)
    {
        std::string ss;
        is >> ss;
        if (ss == "DFS") {
            s = TraverseStrategy::depth_first;
        } else if (ss == "BFS") {
            s = TraverseStrategy::breadth_first;
        } else if (ss == "BFS-VAL") {
            s = TraverseStrategy::breadth_first_value;
        } else if (ss == "UB") {
            s = TraverseStrategy::upper_bound;
        } else {
            throw std::runtime_error("UNKNOWN TRAVERSE STRATEGY");
        }
        return is;
    }


    inline std::istream& operator>>(std::istream& is, BoundStrategy& s)
    {
        std::string ss;
        is >> ss;
        if (ss == "bruteforce") {
            s = BoundStrategy::bruteforce;
        } else if (ss == "local_grob") {
            s = BoundStrategy::local_dual_bound;
        } else if (ss == "local_combined") {
            s = BoundStrategy::local_combined;
        } else if (ss == "local_refined") {
            s = BoundStrategy::local_dual_bound_refined;
        } else if (ss == "local_for_each_point") {
            s = BoundStrategy::local_dual_bound_for_each_point;
        } else {
            throw std::runtime_error("UNKNOWN BOUND STRATEGY");
        }
        return is;
    }

    inline BoundStrategy bs_from_string(std::string s)
    {
        std::stringstream ss(s);
        BoundStrategy result;
        ss >> result;
        return result;
    }

    inline TraverseStrategy ts_from_string(std::string s)
    {
        std::stringstream ss(s);
        TraverseStrategy result;
        ss >> result;
        return result;
    }

    template<class Real>
    struct CalculationParams {
        static constexpr int ALL_DIMENSIONS = -1;

        Real hera_epsilon {0.001}; // relative error in hera call
        Real delta {0.1}; // relative error for matching distance
        int max_depth {8}; // maximal number of refinenemnts
        int initialization_depth {2};
        int dim {0}; // in which dim to calculate the distance; use ALL_DIMENSIONS to get max over all dims
        BoundStrategy bound_strategy {BoundStrategy::local_combined};
        TraverseStrategy traverse_strategy {TraverseStrategy::breadth_first};
        bool tolerate_max_iter_exceeded {false};
        Real actual_error {std::numeric_limits<Real>::max()};
        int actual_max_depth {0};
        int n_hera_calls {0};  // for experiments only; is set in matching_distance function, input value is ignored

        // stop looping over points immediately, if current point's displacement is too large
        // to prune the cell
        // if true, cells are pruned immediately, and bounds may increase
        // (just return something large enough to not prune the cell)
        bool stop_asap { true };

        // print statistics on each quad-tree level
        bool print_stats { false };

#ifdef MD_PRINT_HEAT_MAP
        HeatMaps<Real> heat_maps;
#endif
    };


    template<class Real_, class DiagramProvider>
    class DistanceCalculator {

        using Real = Real_;
        using CellValueVector = std::vector<CellWithValue<Real>>;

    public:
        DistanceCalculator(const DiagramProvider& a,
                const DiagramProvider& b,
                CalculationParams<Real>& params);

        Real distance();

        int get_hera_calls_number() const;

#ifndef MD_TEST_CODE
    private:
#endif

        DiagramProvider module_a_;
        DiagramProvider module_b_;

        CalculationParams<Real>& params_;

        int n_hera_calls_;
        std::map<int, int> n_hera_calls_per_level_;
        Real distance_;

        // if calculate_on_intermediate, then weighted distance
        // will be calculated on centers of each grid in between
        CellValueVector get_refined_grid(int init_depth, bool calculate_on_intermediate, bool calculate_on_last = true);

        CellValueVector get_initial_dual_grid(Real& lower_bound);

#ifdef MD_PRINT_HEAT_MAP
        void heatmap_in_dimension(int dim, int depth);
#endif

        Real get_max_x(int module) const;

        Real get_max_y(int module) const;

        void set_cell_central_value(CellWithValue<Real>& dual_cell);

        Real get_distance();

        Real get_distance_pq();

        Real get_max_possible_value(const CellWithValue<Real>* first_cell_ptr, int n_cells);

        Real get_upper_bound(const CellWithValue<Real>& dual_cell, Real good_enough_upper_bound) const;

        Real get_single_dgm_bound(const CellWithValue<Real>& dual_cell, ValuePoint vp, int module,
                Real good_enough_value) const;

        // this bound depends only on dual box
        Real get_local_dual_bound(int module, const DualBox<Real>& dual_box) const;

        Real get_local_dual_bound(const DualBox<Real>& dual_box) const;

        // this bound depends only on dual box, is more accurate
        Real get_local_refined_bound(int module, const DualBox<Real>& dual_box) const;

        Real get_local_refined_bound(const DualBox<Real>& dual_box) const;

        Real get_good_enough_upper_bound(Real lower_bound) const;

        Real get_max_displacement_single_point(const CellWithValue<Real>& dual_cell, ValuePoint value_point,
                const Point<Real>& p) const;

        void check_upper_bound(const CellWithValue<Real>& dual_cell) const;

        Real distance_on_line(DualPoint<Real> line);
        Real distance_on_line_const(DualPoint<Real> line) const;

        Real current_error(Real lower_bound, Real upper_bound);
    };

    template<class Real>
    Real matching_distance(const Bifiltration<Real>& bif_a, const Bifiltration<Real>& bif_b,
            CalculationParams<Real>& params);

    template<class Real>
    Real matching_distance(const ModulePresentation<Real>& mod_a, const ModulePresentation<Real>& mod_b,
            CalculationParams<Real>& params);

    // for upper bound experiment
    struct UbExperimentRecord {
        double error;
        double lower_bound;
        double upper_bound;
        CellWithValue<double> cell;
        long long int time;
        long long int n_hera_calls;
    };

    inline std::ostream& operator<<(std::ostream& os, const UbExperimentRecord& r)
    {
        os << r.time << "\t" << r.n_hera_calls << "\t" << r.error << "\t" << r.lower_bound << "\t" << r.upper_bound;
        return os;
    }


    template<class K, class V>
    void print_map(const std::map<K, V>& dic)
    {
        for(const auto kv : dic) {
            std::cout << kv.first << " -> " << kv.second << "\n";
        }
    }

} // namespace md

#include "matching/matching_distance.hpp"
