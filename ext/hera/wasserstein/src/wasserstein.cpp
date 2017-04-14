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

#include <assert.h>
#include <algorithm>
#include <functional>
#include <iterator>

#include "def_debug_ws.h"
#include "wasserstein.h"

#ifdef GAUSS_SEIDEL_AUCTION
#include "auction_runner_gs.h"
#else
#include "auction_runner_jac.h"
#endif

namespace geom_ws {

double wassersteinDistVec(const std::vector<DiagramPoint>& A, 
                          const std::vector<DiagramPoint>& B, 
                          const double q, 
                          const double delta,
                          const double _internal_p,
                          const double _initialEpsilon,
                          const double _epsFactor)
{
    if (q < 1) {
#ifndef FOR_R_TDA
        std::cerr << "Wasserstein distance not defined for q = " << q << ", must be >= 1" << std::endl;
#endif
        throw std::runtime_error("Bad q in Wasserstein " + std::to_string(q));
    }
    if (delta < 0.0) {
#ifndef FOR_R_TDA
        std::cerr << "Relative error  " << delta << ", must be > 0" << std::endl;
#endif
        throw std::runtime_error("Bad delta in Wasserstein " + std::to_string(delta));
    }
    if (_initialEpsilon < 0.0) {
#ifndef FOR_R_TDA
        std::cerr << "Initial epsilon = " << _initialEpsilon << ", must be non-negative" << std::endl;
#endif
        throw std::runtime_error("Bad initial epsilon in Wasserstein" + std::to_string(_initialEpsilon));
    }
    if (_epsFactor < 0.0) {
#ifndef FOR_R_TDA
        std::cerr << "Epsilon factor = " << _epsFactor << ", must be non-negative" << std::endl;
#endif
        throw std::runtime_error("Bad epsilon factor in Wasserstein " + std::to_string(_epsFactor));
    }


#ifdef GAUSS_SEIDEL_AUCTION
    AuctionRunnerGS auction(A, B, q,  delta, _internal_p, _initialEpsilon, _epsFactor);
#else
    AuctionRunnerJac auction(A, B, q,  delta, _internal_p);
#endif
    double result = auction.getWassersteinDistance();
    return result;
}

double wassersteinCostVec(const std::vector<DiagramPoint>& A, 
                          const std::vector<DiagramPoint>& B, 
                          const double q, 
                          const double delta,
                          const double _internal_p,
                          const double _initialEpsilon,
                          const double _epsFactor)
{
    if (q < 1) {
#ifndef FOR_R_TDA
        std::cerr << "Wasserstein distance not defined for q = " << q << ", must be >= 1" << std::endl;
#endif
        throw std::runtime_error("Bad q in Wasserstein " + std::to_string(q));
    }
    if (delta < 0.0) {
#ifndef FOR_R_TDA
        std::cerr << "Relative error  " << delta << ", must be > 0" << std::endl;
#endif
        throw std::runtime_error("Bad delta in Wasserstein " + std::to_string(delta));
    }
    if (_initialEpsilon < 0.0) {
#ifndef FOR_R_TDA
        std::cerr << "Initial epsilon = " << _initialEpsilon << ", must be non-negative" << std::endl;
#endif
        throw std::runtime_error("Bad initial epsilon in Wasserstein" + std::to_string(_initialEpsilon));
    }
    if (_epsFactor < 0.0) {
#ifndef FOR_R_TDA
        std::cerr << "Epsilon factor = " << _epsFactor << ", must be non-negative" << std::endl;
#endif
        throw std::runtime_error("Bad epsilon factor in Wasserstein " + std::to_string(_epsFactor));
    }
#ifdef GAUSS_SEIDEL_AUCTION
    AuctionRunnerGS auction(A, B, q,  delta, _internal_p, _initialEpsilon, _epsFactor);
#else
    AuctionRunnerJac auction(A, B, q,  delta, _internal_p);
#endif
    double result = auction.getWassersteinCost();
    return result;
}

bool readDiagramPointSet(const std::string& fname, PairVector& result)
{
    return readDiagramPointSet(fname.c_str(), result);
}

bool readDiagramPointSet(const char* fname, PairVector& result)
{
    size_t lineNumber { 0 };
    result.clear();
    std::ifstream f(fname);
    if (!f.good()) {
#ifndef FOR_R_TDA
        std::cerr << "Cannot open file " << fname << std::endl;
#endif
        return false;
    }
    std::string line;
    while(std::getline(f, line)) {
        lineNumber++;
        // process comments: remove everything after hash
        auto hashPos = line.find_first_of("#", 0);
        if( std::string::npos != hashPos) {
            line = std::string(line.begin(), line.begin() + hashPos);
        }
        if (line.empty()) {
            continue;
        }
         // trim whitespaces 
        auto whiteSpaceFront = std::find_if_not(line.begin(),line.end(),isspace);
        auto whiteSpaceBack = std::find_if_not(line.rbegin(),line.rend(),isspace).base();
        if (whiteSpaceBack <= whiteSpaceFront) {
            // line consists of spaces only - move to the next line
            continue;
        }
        line = std::string(whiteSpaceFront,whiteSpaceBack);
        double x, y;
        std::istringstream iss(line);
        if (not(iss >> x >> y)) {
#ifndef FOR_R_TDA
            std::cerr << "Error in file " << fname << ", line number " << lineNumber << ": cannot parse \"" << line << "\"" << std::endl;
#endif
            return false;
        }
        result.push_back(std::make_pair(x,y));
    }
    f.close();
    return true;
}


void removeDuplicates(PairVector& dgmA, PairVector& dgmB)
{
    std::map<std::pair<double, double>, int> mapA, mapB;
    // copy points to maps
    for(const auto& ptA : dgmA) {
        mapA[ptA]++;
    }
    for(const auto& ptB : dgmB) {
        mapB[ptB]++;
    }
    // clear vectors
    dgmA.clear();
    dgmB.clear();
    // remove duplicates from maps
    for(auto& pointMultiplicityPair : mapA) {
        auto iterB = mapB.find(pointMultiplicityPair.first);
        if (iterB != mapB.end()) {
            int duplicateMultiplicity = std::min(pointMultiplicityPair.second, iterB->second);
            pointMultiplicityPair.second -= duplicateMultiplicity;
            iterB->second -= duplicateMultiplicity;
        }
    }
    // copy points back to vectors
    for(const auto& pointMultiplicityPairA : mapA) {
        assert( pointMultiplicityPairA.second >= 0);
        for(int i = 0; i < pointMultiplicityPairA.second; ++i) {
            dgmA.push_back(pointMultiplicityPairA.first);
        }
    }

    for(const auto& pointMultiplicityPairB : mapB) {
        assert( pointMultiplicityPairB.second >= 0);
        for(int i = 0; i < pointMultiplicityPairB.second; ++i) {
            dgmB.push_back(pointMultiplicityPairB.first);
        }
    }
}

} // end of namespace geom_ws
