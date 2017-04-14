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

#ifndef WASSERSTEIN_H
#define WASSERSTEIN_H

#include <vector>
#include <map>
#include <math.h>

#include "basic_defs_ws.h"

// use Gauss-Seidel version; comment out to switch to Jacobi (not recommended)
#define GAUSS_SEIDEL_AUCTION

namespace geom_ws {

using PairVector = std::vector<std::pair<double, double>>;

// get Wasserstein distance between two persistence diagrams
double wassersteinDistVec(const std::vector<DiagramPoint>& A, 
                          const std::vector<DiagramPoint>& B, 
                          const double q, 
                          const double delta,
                          const double _internal_p = std::numeric_limits<double>::infinity(),
                          const double _initialEpsilon = 0.0,
                          const double _epsFactor = 0.0);

// get Wasserstein cost (distance^q) between two persistence diagrams
double wassersteinCostVec(const std::vector<DiagramPoint>& A, 
                          const std::vector<DiagramPoint>& B, 
                          const double q, 
                          const double delta,
                          const double _internal_p = std::numeric_limits<double>::infinity(),
                          const double _initialEpsilon = 0.0,
                          const double _epsFactor = 0.0);


// compare as multisets
template<class PairContainer>
bool areEqual(PairContainer& dgm1, PairContainer& dgm2)
{
    if (dgm1.size() != dgm2.size()) {
        return false;
    }

    std::map<std::pair<double, double>, int> m1, m2;

    for(const auto& pair1 : dgm1) {
        m1[pair1]++;
    }

    for(const auto& pair2 : dgm2) {
        m2[pair2]++;
    }

    return m1 == m2;
}

template<class PairContainer>
double wassersteinDist(PairContainer& A, PairContainer& B, const double q, const double delta, const double _internal_p = std::numeric_limits<double>::infinity(), const double _initialEpsilon = 0.0, const double _epsFactor = 0.0)
{
    if (areEqual(A, B)) {
        return 0.0;
    }

    std::vector<DiagramPoint> dgmA, dgmB;
    // loop over A, add projections of A-points to corresponding positions
    // in B-vector
    for(auto& pairA : A) {
        double x = pairA.first;
        double y = pairA.second;
        dgmA.push_back(DiagramPoint(x, y,  DiagramPoint::NORMAL));
        dgmB.push_back(DiagramPoint(x, y,  DiagramPoint::DIAG));
    }
    // the same for B
    for(auto& pairB : B) {
        double x = pairB.first;
        double y = pairB.second;
        dgmA.push_back(DiagramPoint(x, y,  DiagramPoint::DIAG));
        dgmB.push_back(DiagramPoint(x, y,  DiagramPoint::NORMAL));
    }
    
    return wassersteinDistVec(dgmA, dgmB, q, delta, _internal_p, _initialEpsilon, _epsFactor);
}

template<class PairContainer>
double wassersteinCost(PairContainer& A, PairContainer& B, const double q, const double delta, const double _internal_p = std::numeric_limits<double>::infinity(), const double _initialEpsilon = 0.0, const double _epsFactor = 0.0)
{
    if (areEqual(A, B)) {
        return 0.0;
    }

    std::vector<DiagramPoint> dgmA, dgmB;
    // loop over A, add projections of A-points to corresponding positions
    // in B-vector
    for(auto& pairA : A) {
        double x = pairA.first;
        double y = pairA.second;
        dgmA.push_back(DiagramPoint(x, y,  DiagramPoint::NORMAL));
        dgmB.push_back(DiagramPoint(x, y,  DiagramPoint::DIAG));
    }
    // the same for B
    for(auto& pairB : B) {
        double x = pairB.first;
        double y = pairB.second;
        dgmA.push_back(DiagramPoint(x, y,  DiagramPoint::DIAG));
        dgmB.push_back(DiagramPoint(x, y,  DiagramPoint::NORMAL));
    }
    
    return wassersteinCostVec(dgmA, dgmB, q, delta, _internal_p, _initialEpsilon, _epsFactor);
}


// fill in result with points from file fname
// return false if file can't be opened
// or error occurred while reading
bool readDiagramPointSet(const char* fname, PairVector& result);
bool readDiagramPointSet(const std::string& fname, PairVector& result);
void removeDuplicates(PairVector& dgmA, PairVector& dgmB);
 
} // end of namespace geom_ws

#endif
