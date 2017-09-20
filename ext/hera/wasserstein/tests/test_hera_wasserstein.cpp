#define LOG_AUCTION
#include "catch/catch.hpp"

#include <sstream>
#include <iostream>


#include "wasserstein.h"



using PairVector = std::vector<std::pair<double, double>>;

std::vector<std::string> split_on_delim(const std::string& s, char delim)
{
    std::stringstream ss(s);
    std::string token;
    std::vector<std::string> tokens;
    while(std::getline(ss, token, delim)) {
        tokens.push_back(token);
    }
    return tokens;
}


// single row in a file with test cases
struct TestFromFileCase {

    std::string file_1;
    std::string file_2;
    double q;
    double delta;
    double internal_p;

    TestFromFileCase(std::string s)
{
    auto tokens = split_on_delim(s, ' ');
    assert(tokens.size() == 5);

    file_1 = tokens[0];
    file_2 = tokens[1];
    q = std::stod(tokens[2]);
    delta = std::stod(tokens[3]);
    internal_p = std::stod(tokens[4]);

    if ( q < 1.0 or std::isinf(q) or delta <= 0.0 or
        (internal_p != hera::get_infinity<double>() and internal_p < 1.0)) {
        throw std::runtime_error("Bad line in test_list.txt");
    }
}
};


TEST_CASE("simple cases", "wasserstein_dist")
{
    PairVector diagram_A, diagram_B;
    hera::AuctionParams<double> params;
    params.wasserstein_power = 1.0;
    params.delta = 0.01;
    params.internal_p = hera::get_infinity<double>();
    params.initial_epsilon = 0.0;
    params.epsilon_common_ratio = 0.0;
    params.max_num_phases = 30;
    params.gamma_threshold = 0.0;
    params.max_bids_per_round = 0;  // use Jacobi

    SECTION("trivial: two empty diagrams") {
        REQUIRE(  0.0 == hera::wasserstein_dist<>(diagram_A, diagram_B, params));
    }

    SECTION("trivial: one empty diagram, one single-point diagram") {

        diagram_A.emplace_back(1.0, 2.0);

        double d1 = hera::wasserstein_dist<>(diagram_A, diagram_B, params);
        REQUIRE(  fabs(d1 - 0.5) <= 0.00000000001 );

        double d2 = hera::wasserstein_dist<>(diagram_B, diagram_A, params);
        REQUIRE(  fabs(d2 - 0.5) <= 0.00000000001 );

        params.internal_p = 2.0;
        double corr_answer = 1.0 / std::sqrt(2.0);
        double d3 = hera::wasserstein_dist<>(diagram_B, diagram_A, params);
        REQUIRE(  fabs(d3 - corr_answer) <= 0.00000000001 );

    }

    SECTION("trivial: two single-point diagrams-1") {

        diagram_A.emplace_back(10.0, 20.0);  // (5, 5)
        diagram_B.emplace_back(13.0, 19.0);  // (3, 3)

        double d1 = hera::wasserstein_dist<>(diagram_A, diagram_B, params);
        double d2 = hera::wasserstein_dist<>(diagram_B, diagram_A, params);
        REQUIRE(  fabs(d1 - d2) <= 0.00000000001 );
        REQUIRE(  fabs(d1 - 3.0) <= 0.00000000001 );

        params.wasserstein_power = 2.0;
        double d3 = hera::wasserstein_cost<>(diagram_A, diagram_B, params);
        double d4 = hera::wasserstein_cost<>(diagram_B, diagram_A, params);
        REQUIRE(  fabs(d3 - d4) <= 0.00000000001 );
        REQUIRE(  fabs(d4 - 9.0) <= 0.00000000001 );

        params.wasserstein_power = 1.0;
        params.internal_p = 1.0;
        double d5 = hera::wasserstein_cost<>(diagram_A, diagram_B, params);
        double d6 = hera::wasserstein_cost<>(diagram_B, diagram_A, params);
        REQUIRE(  fabs(d5 - d6) <= 0.00000000001 );
        REQUIRE(  fabs(d5 - 4.0) <= 0.00000000001 );

        params.internal_p = 2.0;
        double d7 = hera::wasserstein_cost<>(diagram_A, diagram_B, params);
        double d8 = hera::wasserstein_cost<>(diagram_B, diagram_A, params);
        REQUIRE(  fabs(d7 - d8) <= 0.00000000001 );
        REQUIRE(  fabs(d7 - std::sqrt(10.0)) <= 0.00000000001 );

    }

    SECTION("trivial: two single-point diagrams-2") {

        diagram_A.emplace_back(10.0, 20.0);  // (5, 5)
        diagram_B.emplace_back(130.0, 138.0);  // (4, 4)

        double d1 = hera::wasserstein_cost<>(diagram_A, diagram_B, params);
        double d2 = hera::wasserstein_cost<>(diagram_B, diagram_A, params);
        REQUIRE(  fabs(d1 - d2) <= 0.00000000001 );
        REQUIRE(  fabs(d1 - 9.0) <= 0.00000000001 );

        params.wasserstein_power = 2.0;
        double d3 = hera::wasserstein_cost<>(diagram_A, diagram_B, params);
        double d4 = hera::wasserstein_cost<>(diagram_B, diagram_A, params);
        REQUIRE(  fabs(d3 - d4) <= 0.00000000001 );
        REQUIRE(  fabs(d4 - 41.0) <= 0.00000000001 ); // 5^2 + 4^2

        params.wasserstein_power = 1.0;
        params.internal_p = 1.0;
        double d5 = hera::wasserstein_cost<>(diagram_A, diagram_B, params);
        double d6 = hera::wasserstein_cost<>(diagram_B, diagram_A, params);
        REQUIRE(  fabs(d5 - d6) <= 0.00000000001 );
        REQUIRE(  fabs(d5 - 18.0) <= 0.00000000001 ); // 5 + 5 + 4 + 4

        params.internal_p = 2.0;
        double d7 = hera::wasserstein_cost<>(diagram_A, diagram_B, params);
        double d8 = hera::wasserstein_cost<>(diagram_B, diagram_A, params);
        REQUIRE(  fabs(d7 - d8) <= 0.00000000001 );
        REQUIRE(  fabs(d7 - 9 * std::sqrt(2.0)) <= 0.00000000001 ); // sqrt(5^2 + 5^2) + sqrt(4^2 + 4^2) = 9 sqrt(2)

    }



    SECTION("from file:") {
        const char* file_name = "../tests/data/test_list.txt";
        std::ifstream f;
        f.open(file_name);
        std::vector<TestFromFileCase> test_params;
        std::string s;
        while (std::getline(f, s)) {
            test_params.emplace_back(s);
        }
    }



}
