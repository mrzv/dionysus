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

#ifndef HERA_BOTTLENECK_HPP
#define HERA_BOTTLENECK_HPP

#ifdef FOR_R_TDA
#undef DEBUG_BOUND_MATCH
#undef DEBUG_MATCHING
#undef VERBOSE_BOTTLENECK
#endif


#include <iomanip>
#include <sstream>
#include <string>
#include <cctype>

#include "bottleneck_detail.h"

namespace hera {
namespace bt {

// return the interval (distMin, distMax) such that:
// a) actual bottleneck distance between A and B is contained in the interval
// b) if the interval is not (0,0), then  (distMax - distMin) / distMin < epsilon
template<class Real>
std::pair<Real, Real> bottleneckDistApproxInterval(DiagramPointSet<Real>& A, DiagramPointSet<Real>& B, const Real epsilon)
{
    // empty diagrams are not considered as error
    if (A.empty() and B.empty())
        return std::make_pair(0.0, 0.0);

    // link diagrams A and B by adding projections
    addProjections(A, B);

    // TODO: think about that!
    // we need one threshold for checking if the distance is 0,
    // another one for the oracle!
    constexpr Real epsThreshold { 1.0e-10 };
    std::pair<Real, Real> result { 0.0, 0.0 };
    bool useRangeSearch { true };
    // construct an oracle
    BoundMatchOracle<Real> oracle(A, B, epsThreshold, useRangeSearch);
    // check for distance = 0
    if (oracle.isMatchLess(2*epsThreshold)) {
        return result;
    }
    // get a 3-approximation of maximal distance between A and B
    // as a starting value for probe distance
    Real distProbe  { getFurthestDistance3Approx<Real, DiagramPointSet<Real>>(A, B) };
    // aliases for result components
    Real& distMin {result.first};
    Real& distMax {result.second};

    if ( oracle.isMatchLess(distProbe) ) {
        // distProbe is an upper bound,
        // find lower bound with binary search
        do {
            distMax = distProbe;
            distProbe /= 2.0;
        } while (oracle.isMatchLess(distProbe));
        distMin = distProbe;
    } else {
        // distProbe is a lower bound,
        // find upper bound with exponential search
        do {
            distMin = distProbe;
            distProbe *= 2.0;
        } while (!oracle.isMatchLess(distProbe));
        distMax = distProbe;
    }
    // bounds are found, perform binary search
    //std::cout << "Bounds found, distMin = " << distMin << ", distMax = " << distMax << ", ratio = " << ( distMax - distMin ) / distMin << std::endl ;
    distProbe = ( distMin + distMax ) / 2.0;
    while (  ( distMax - distMin ) / distMin >= epsilon ) {
        if (oracle.isMatchLess(distProbe)) {
            distMax = distProbe;
        } else {
            distMin = distProbe;
        }
        distProbe = ( distMin + distMax ) / 2.0;
    }
    return result;
}

template<class Real>
void sampleDiagramForHeur(const DiagramPointSet<Real>& dgmIn, DiagramPointSet<Real>& dgmOut)
{
#ifdef VERBOSE_BOTTLENECK
    std::cout << "Entered sampleDiagramForHeur, dgmIn.size = " << dgmIn.size() << std::endl;
#endif
    struct pair_hash {
        std::size_t operator()(const std::pair<Real, Real> p) const
        {
            return std::hash<Real>()(p.first) ^ std::hash<Real>()(p.second);
        }
    };
    std::unordered_map<std::pair<Real, Real>, int, pair_hash> m;
    for(auto ptIter = dgmIn.cbegin(); ptIter != dgmIn.cend(); ++ptIter) {
        if (ptIter->isNormal()) {
            m[std::make_pair(ptIter->getRealX(), ptIter->getRealY())]++;
        }
    }
#ifdef VERBOSE_BOTTLENECK
    std::cout << "map filled in, m.size = " << m.size() << std::endl;
#endif
     if (m.size() < 2) {
        dgmOut = dgmIn;
        return;
    }
    std::vector<int> v;
    for(const auto& ptQtyPair : m) {
        v.push_back(ptQtyPair.second);
    }
#ifdef VERBOSE_BOTTLENECK
    std::cout << "v filled in, v.size = " << v.size() << std::endl;
#endif
    std::sort(v.begin(), v.end());
#ifdef VERBOSE_BOTTLENECK
    std::cout << "v sorted" << std::endl;
#endif
    int maxLeap = v[1] - v[0];
    int cutVal = v[0];
    for(int i = 1; i < static_cast<int>(v.size())- 1; ++i) {
        int currLeap = v[i+1] - v[i];
        if (currLeap > maxLeap) {
            maxLeap = currLeap;
            cutVal = v[i];
        }
    }
#ifdef VERBOSE_BOTTLENECK
    std::cout << "cutVal found, cutVal = " << cutVal << std::endl;
#endif
     std::vector<std::pair<Real, Real>> vv;
    // keep points whose multiplicites are at most cutVal
    // quick-and-dirty: fill in vv with copies of each point
    // to construct DiagramPointSet from it later
    for(const auto& ptQty : m) {
        if (ptQty.second < cutVal) {
            for(int i = 0; i < ptQty.second; ++i) {
                vv.push_back(std::make_pair(ptQty.first.first, ptQty.first.second));
            }
        }
    }
#ifdef VERBOSE_BOTTLENECK
    std::cout << "vv filled in, vv.size = " << v.size() << std::endl;
#endif
    dgmOut.clear();
    dgmOut = DiagramPointSet<Real>(vv.begin(), vv.end());
#ifdef VERBOSE_BOTTLENECK
    std::cout << "dgmOut filled in, dgmOut.size = " << dgmOut.size() << std::endl;
#endif
}


// return the interval (distMin, distMax) such that:
// a) actual bottleneck distance between A and B is contained in the interval
// b) if the interval is not (0,0), then  (distMax - distMin) / distMin < epsilon
template<class Real>
std::pair<Real, Real> bottleneckDistApproxIntervalWithInitial(DiagramPointSet<Real>& A, DiagramPointSet<Real>& B, const Real epsilon, const std::pair<Real, Real> initialGuess)
{
    // empty diagrams are not considered as error
    if (A.empty() and B.empty())
        return std::make_pair(0.0, 0.0);

    // link diagrams A and B by adding projections
    addProjections(A, B);

    constexpr Real epsThreshold { 1.0e-10 };
    std::pair<Real, Real> result { 0.0, 0.0 };
    bool useRangeSearch { true };
    // construct an oracle
    BoundMatchOracle<Real> oracle(A, B, epsThreshold, useRangeSearch);
    Real& distMin {result.first};
    Real& distMax {result.second};

    // initialize search interval from initialGuess
    distMin = initialGuess.first;
    distMax = initialGuess.second;

    assert(distMin <= distMax);

    // make sure that distMin is a lower bound
    while(oracle.isMatchLess(distMin)) {
        // distMin is in fact an upper bound, so assign it to distMax
        distMax = distMin;
        // and decrease distMin by 5 %
        distMin = 0.95 * distMin;
    }

    // make sure that distMax is an upper bound
    while(not oracle.isMatchLess(distMax)) {
        // distMax is in fact a lower bound, so assign it to distMin
        distMin = distMax;
        // and increase distMax by 5 %
        distMax = 1.05 * distMax;
    }

    // bounds are found, perform binary search
    //std::cout << "Bounds found, distMin = " << distMin << ", distMax = " << distMax << ", ratio = " << ( distMax - distMin ) / distMin << std::endl ;
    Real distProbe = ( distMin + distMax ) / 2.0;
    while (  ( distMax - distMin ) / distMin >= epsilon ) {
        if (oracle.isMatchLess(distProbe)) {
            distMax = distProbe;
        } else {
            distMin = distProbe;
        }
        distProbe = ( distMin + distMax ) / 2.0;
    }
    return result;
}

// return the interval (distMin, distMax) such that:
// a) actual bottleneck distance between A and B is contained in the interval
// b) if the interval is not (0,0), then  (distMax - distMin) / distMin < epsilon
// use heuristic: initial estimate on sampled diagrams
template<class Real>
std::pair<Real, Real> bottleneckDistApproxIntervalHeur(DiagramPointSet<Real>& A, DiagramPointSet<Real>& B, const Real epsilon)
{
    // empty diagrams are not considered as error
    if (A.empty() and B.empty())
        return std::make_pair(0.0, 0.0);

    DiagramPointSet<Real> sampledA, sampledB;
    sampleDiagramForHeur(A, sampledA);
    sampleDiagramForHeur(B, sampledB);
#ifdef VERBOSE_BOTTLENECK
    std::cout << "A : " << A.size() << ", sampled: " << sampledA.size() << std::endl;
    std::cout << "B : " << B.size() << ", sampled: " << sampledB.size() << std::endl;
#endif
    std::pair<Real, Real> initGuess = bottleneckDistApproxInterval(sampledA, sampledB, epsilon);
#ifdef VERBOSE_BOTTLENECK
    std::cout << "initial guess with sampling: " << initGuess.first << ", " << initGuess.second << std::endl;
    std::cout << "running on the original diagrams" << std::endl;
#endif
    return bottleneckDistApproxIntervalWithInitial<Real>(A, B, epsilon, initGuess);
}



// get approximate distance,
// see bottleneckDistApproxInterval
template<class Real>
Real bottleneckDistApprox(DiagramPointSet<Real>& A, DiagramPointSet<Real>& B, const Real epsilon)
{
    auto interval = bottleneckDistApproxInterval<Real>(A, B, epsilon);
    return interval.second;
}


template<class Real>
Real bottleneckDistExactFromSortedPwDist(DiagramPointSet<Real>&A, DiagramPointSet<Real>& B, std::vector<Real>& pairwiseDist, const int decPrecision)
{
    //for(size_t k = 0; k < pairwiseDist.size(); ++k) {
        //std::cout << "pairwiseDist[" << k << "] = " << std::setprecision(15) << pairwiseDist[k] << std::endl;
    //}
    // trivial case: we have only one candidate
    if (pairwiseDist.size() == 1)
        return pairwiseDist[0];

    bool useRangeSearch = true;
    Real distEpsilon = std::numeric_limits<Real>::max();
    Real diffThreshold = 0.1;
    for(int k = 0; k < decPrecision; ++k) {
        diffThreshold /= 10.0;
    }
    for(size_t k = 0; k < pairwiseDist.size() - 2; ++k) {
        auto diff = pairwiseDist[k+1]- pairwiseDist[k];
        //std::cout << "diff = " << diff << ", pairwiseDist[k] = " << pairwiseDist[k] << std::endl;
        if ( diff > diffThreshold and diff < distEpsilon ) {
            distEpsilon = diff;
        }
    }
    distEpsilon /= 3.0;
    //std::cout << "decPrecision = " << decPrecision << ", distEpsilon = " << distEpsilon << std::endl;

    BoundMatchOracle<Real> oracle(A, B, distEpsilon, useRangeSearch);
    // binary search
    size_t iterNum {0};
    size_t idxMin {0}, idxMax {pairwiseDist.size() - 1};
    size_t idxMid;
    while(idxMax > idxMin) {
        idxMid = static_cast<size_t>(floor(idxMin + idxMax) / 2.0);
        //std::cout << "while begin: min = " << idxMin << ", idxMax = " << idxMax << ", idxMid = " << idxMid << ", testing d = " << std::setprecision(15) << pairwiseDist[idxMid] << std::endl;
        iterNum++;
        // not A[imid] < dist <=>  A[imid] >= dist  <=> A[imid[ >= dist + eps
        if (oracle.isMatchLess(pairwiseDist[idxMid] + distEpsilon / 2.0)) {
            //std::cout << "isMatchLess = true" << std::endl;
            idxMax = idxMid;
        } else {
            //std::cout << "isMatchLess = false " << std::endl;
            idxMin = idxMid + 1;
        }
        //std::cout << "while end: idxMin = " << idxMin << ", idxMax = " << idxMax << ", idxMid = " << idxMid << std::endl;
    }
    idxMid = static_cast<size_t>(floor(idxMin + idxMax) / 2.0);
    return pairwiseDist[idxMid];
}


template<class Real>
Real bottleneckDistExact(DiagramPointSet<Real>& A, DiagramPointSet<Real>& B)
{
    return bottleneckDistExact(A, B, 14);
}

template<class Real>
Real bottleneckDistExact(DiagramPointSet<Real>& A, DiagramPointSet<Real>& B, const int decPrecision)
{
    using DgmPoint = DiagramPoint<Real>;

    constexpr Real epsilon = 0.001;
    auto interval = bottleneckDistApproxInterval(A, B, epsilon);
    const Real delta = 0.50001 * (interval.second - interval.first);
    const Real approxDist = 0.5 * ( interval.first + interval.second);
    const Real minDist = interval.first;
    const Real maxDist = interval.second;
    //std::cout << std::setprecision(15) <<  "minDist = " << minDist << ", maxDist = " << maxDist << std::endl;
    if ( delta == 0 ) {
        return interval.first;
    }
    // copy points from A to a vector
    // todo: get rid of this?
    std::vector<DgmPoint> pointsA;
    pointsA.reserve(A.size());
    for(const auto& ptA : A) {
        pointsA.push_back(ptA);
    }

    //std::vector<Real> killdist;
    //for(auto pta : a) {
        //for(auto ptb : b) {
            //if ( distlinf(pta, ptb) > mindist and distlinf(pta, ptb) < maxdist) {
                //killdist.push_back(distlinf(pta, ptb));
                //std::cout << pta << ", " << ptb << std::endl;
            //}
        //}
    //}
    //std::sort(killdist.begin(), killdist.end());
    //for(auto d : killdist) {
        //std::cout << d << std::endl;
    //}
    //std::cout << "*************" << std::endl;

    // in this vector we store the distances between the points
    // that are candidates to realize
    std::vector<Real> pairwiseDist;
    {
        // vector to store centers of vertical stripes
        // two for each point in A and the id of the corresponding point
        std::vector<std::pair<Real, DgmPoint>> xCentersVec;
        xCentersVec.reserve(2 * pointsA.size());
        for(auto ptA : pointsA) {
            xCentersVec.push_back(std::make_pair(ptA.getRealX() - approxDist, ptA));
            xCentersVec.push_back(std::make_pair(ptA.getRealX() + approxDist, ptA));
        }
        // lambda to compare pairs <coordinate, id> w.r.t coordinate
        auto compLambda = [](std::pair<Real, DgmPoint> a, std::pair<Real, DgmPoint> b)
                                { return a.first < b.first; };

        std::sort(xCentersVec.begin(), xCentersVec.end(), compLambda);
        //std::cout << "xCentersVec.size  = " << xCentersVec.size() << std::endl;
        //for(auto p = xCentersVec.begin(); p!= xCentersVec.end(); ++p) {
            //if (p->second.id == 200) {
                //std::cout << "index of 200: " << p - xCentersVec.begin() << std::endl;
            //}
        //}
        //std::vector<DgmPoint>
        // todo: sort points in B, reduce search range in lower and upper bounds
        for(auto ptB : B) {
            // iterator to the first stripe such that ptB lies to the left
            // from its right boundary (x_B <= x_j + \delta iff x_j >= x_B - \delta
            auto itStart = std::lower_bound(xCentersVec.begin(),
                                            xCentersVec.end(),
                                            std::make_pair(ptB.getRealX() - delta, ptB),
                                            compLambda);
            //if (ptB.id == 236) {
                //std::cout << itStart - xCentersVec.begin() <<  std::endl;
            //}

            for(auto iterA = itStart; iterA < xCentersVec.end(); ++iterA) {
                //if (ptB.id == 236) {
                    //std::cout << "consider " << iterA->second << std::endl;
                //}
                if ( ptB.getRealX() < iterA->first - delta) {
                    // from that moment x_B >= x_j - delta
                    // is violated: x_B no longer lies to right from the left
                    // boundary of current stripe
                    //if (ptB.id == 236) {
                        //std::cout << "break" << std::endl;
                    //}
                    break;
                }
                // we're here => ptB lies in vertical stripe,
                // check if distance fits into the interval we've found
                Real pwDist = distLInf(iterA->second, ptB);
                //if (ptB.id == 236) {
                    //std::cout << pwDist << std::endl;
                //}
                //std::cout << 1000*minDist << " <= " << 1000*pwDist << " <= " << 1000*maxDist << std::endl;
                if (pwDist >= minDist and pwDist <= maxDist) {
                    pairwiseDist.push_back(pwDist);
                }
            }
        }
    }

    {
        // for y
        // vector to store centers of vertical stripes
        // two for each point in A and the id of the corresponding point
        std::vector<std::pair<Real, DgmPoint>> yCentersVec;
        yCentersVec.reserve(2 * pointsA.size());
        for(auto ptA : pointsA) {
            yCentersVec.push_back(std::make_pair(ptA.getRealY() - approxDist, ptA));
            yCentersVec.push_back(std::make_pair(ptA.getRealY() + approxDist, ptA));
        }
        // lambda to compare pairs <coordinate, id> w.r.t coordinate
        auto compLambda = [](std::pair<Real, DgmPoint> a, std::pair<Real, DgmPoint> b)
                                { return a.first < b.first; };

        std::sort(yCentersVec.begin(), yCentersVec.end(), compLambda);

        // std::cout << "Sorted vector of y-centers:" << std::endl;
        //for(auto coordPtPair : yCentersVec) {
            //std::cout << coordPtPair.first << ", id = " << coordPtPair.second.id << std::endl;
        //}
        /*std::cout << "End of sorted vector of y-centers:" << std::endl;*/

        //std::vector<DgmPoint>
        // todo: sort points in B, reduce search range in lower and upper bounds
        for(auto ptB : B) {
            auto itStart = std::lower_bound(yCentersVec.begin(),
                                            yCentersVec.end(),
                                            std::make_pair(ptB.getRealY() - delta, ptB),
                                            compLambda);


            for(auto iterA = itStart; iterA < yCentersVec.end(); ++iterA) {
                if ( ptB.getRealY() < iterA->first - delta) {
                    break;
                }
                Real pwDist = distLInf(iterA->second, ptB);
                //std::cout << 1000*minDist << " <= " << 1000*pwDist << " <= " << 1000*maxDist << std::endl;
                if (pwDist >= minDist and pwDist <= maxDist) {
                    pairwiseDist.push_back(pwDist);
                }
            }
        }
    }

    //std::cout << "pairwiseDist.size = " << pairwiseDist.size() << " out of " << A.size() * A.size() << std::endl;
    std::sort(pairwiseDist.begin(), pairwiseDist.end());
    //for(auto ddd : pairwiseDist) {
        //std::cout << std::setprecision(15) << ddd << std::endl;
    //}

    return bottleneckDistExactFromSortedPwDist(A, B, pairwiseDist, decPrecision);
}

template<class Real>
Real bottleneckDistSlow(DiagramPointSet<Real>& A, DiagramPointSet<Real>& B)
{
    using DistVerticesPair = std::pair<Real, std::pair<size_t, size_t>>;

    // use range search when building the layer graph
    bool useRangeSearch { true };
    // find maximum of min. distances for each point,
    // use this value as lower bound for bottleneck distance
    bool useHeurMinIdx { true };

    // find matching in a greedy manner to
    // get an upper bound for a bottleneck distance
    bool useHeurGreedyMatching { false };

    // use successive multiplication of idxMin with 2 to get idxMax
    bool goUpToFindIdxMax { false };
    //
    goUpToFindIdxMax = goUpToFindIdxMax and !useHeurGreedyMatching;

    if (!useHeurGreedyMatching) {
        long int N = 3 * (A.size() / 2 ) * (B.size() / 2);
        std::vector<Real> pairwiseDist;
        pairwiseDist.reserve(N);
        Real maxMinDist {0.0};
        for(auto& p_A : A) {
            Real minDist { std::numeric_limits<Real>::max() };
            for(auto& p_B : B) {
                if (p_A.isNormal() or p_B.isNormal()) {
                    Real d = distLInf(p_A, p_B);
                    pairwiseDist.push_back(d);
                    if (useHeurMinIdx and p_A.isNormal()) {
                        if (d < minDist)
                            minDist = d;
                    }
                }
            }
            if (useHeurMinIdx and p_A.isNormal() and minDist > maxMinDist) {
                maxMinDist = minDist;
            }
        }

        std::sort(pairwiseDist.begin(), pairwiseDist.end());

        Real distEpsilon = std::numeric_limits<Real>::max();
        for(size_t k = 0; k < pairwiseDist.size() - 2; ++k) {
            auto diff = pairwiseDist[k+1]- pairwiseDist[k];
            if ( diff > 1.0e-10 and diff < distEpsilon ) {
                distEpsilon = diff;
            }
        }
        distEpsilon /= 3.0;

        BoundMatchOracle<Real> oracle(A, B, distEpsilon, useRangeSearch);
        // binary search
        size_t iterNum {0};
        size_t idxMin {0}, idxMax {pairwiseDist.size() - 1};
        if (useHeurMinIdx) {
            auto maxMinIter = std::equal_range(pairwiseDist.begin(), pairwiseDist.end(), maxMinDist);
            assert(maxMinIter.first != pairwiseDist.end());
            idxMin = maxMinIter.first - pairwiseDist.begin();
            //std::cout << "maxMinDist = " << maxMinDist << ", idxMin = " << idxMin << ", d = " << pairwiseDist[idxMin] << std::endl;
        }

        if (goUpToFindIdxMax) {
            if ( pairwiseDist.size() == 1) {
                return pairwiseDist[0];
            }

            idxMax = std::max<size_t>(idxMin, 1);
            while (!oracle.isMatchLess(pairwiseDist[idxMax])) {
                //std::cout << "entered while" << std::endl;
                idxMin = idxMax;
                if (2*idxMax > pairwiseDist.size() -1) {
                    idxMax = pairwiseDist.size() - 1;
                    break;
                } else {
                    idxMax *= 2;
                }
            }
            //std::cout << "size = " << pairwiseDist.size() << ", idxMax = " << idxMax <<  ", pw[max] = " << pairwiseDist[idxMax] << std::endl;
        }

        size_t idxMid { (idxMin + idxMax) / 2 };
        while(idxMax > idxMin) {
            iterNum++;
            if (oracle.isMatchLess(pairwiseDist[idxMid])) {
                idxMax = idxMid;
            } else {
                if (idxMax - idxMin == 1)
                    idxMin++;
                else
                    idxMin = idxMid;
            }
            idxMid = (idxMin + idxMax) / 2;
        }
        return pairwiseDist[idxMid];
    } else {
        // with greeedy matching
        long int N = A.size() * B.size();
        std::vector<DistVerticesPair> pairwiseDist;
        pairwiseDist.reserve(N);
        Real maxMinDist {0.0};
        size_t idxA{0}, idxB{0};
        for(auto p_A : A) {
            Real minDist { std::numeric_limits<Real>::max() };
            idxB = 0;
            for(auto p_B : B) {
                Real d = distLInf(p_A, p_B);
                pairwiseDist.push_back( std::make_pair(d, std::make_pair(idxA, idxB) ) );
                if (useHeurMinIdx and p_A.isNormal()) {
                    if (d < minDist)
                        minDist = d;
                }
                idxB++;
            }
            if (useHeurMinIdx and p_A.isNormal() and minDist > maxMinDist) {
                maxMinDist = minDist;
            }
            idxA++;
        }

        auto compLambda = [](DistVerticesPair a, DistVerticesPair b)
                    { return a.first < b.first;};

        std::sort(pairwiseDist.begin(),
                  pairwiseDist.end(),
                  compLambda);

        Real distEpsilon = std::numeric_limits<Real>::max();
        for(size_t k = 0; k < pairwiseDist.size() - 2; ++k) {
            auto diff = pairwiseDist[k+1].first - pairwiseDist[k].first;
            if ( diff > 1.0e-10 and diff < distEpsilon ) {
                distEpsilon = diff;
            }
        }
        distEpsilon /= 3.0;

        BoundMatchOracle<Real> oracle(A, B, distEpsilon, useRangeSearch);

        // construct greedy matching
        size_t numVert { A.size() };
        size_t numMatched { 0 };
        std::unordered_set<size_t> aTobMatched, bToaMatched;
        aTobMatched.reserve(numVert);
        bToaMatched.reserve(numVert);
        size_t distVecIdx {0};
        while( numMatched < numVert) {
            auto vertPair = pairwiseDist[distVecIdx++].second;
            //std::cout << "distVecIdx = " << distVecIdx <<   ", matched: " << numMatched << " out of " << numVert << std::endl;
            //std::cout << "vertex A idx = " << vertPair.first <<   ", B idx: " << vertPair.second << " out of " << numVert << std::endl;
            if ( aTobMatched.count(vertPair.first) == 0 and
                 bToaMatched.count(vertPair.second) == 0 ) {
                aTobMatched.insert(vertPair.first);
                bToaMatched.insert(vertPair.second);
                numMatched++;
            }
        }
        size_t idxMax = distVecIdx-1;
        //std::cout << "idxMax = " << idxMax << ", size = " << pairwiseDist.size() << std::endl;
        // binary search
        size_t iterNum {0};
        size_t idxMin {0};
        if (useHeurMinIdx) {
            auto maxMinIter = std::equal_range(pairwiseDist.begin(),
                                               pairwiseDist.end(),
                                               std::make_pair(maxMinDist, std::make_pair(0,0)),
                                               compLambda);
            assert(maxMinIter.first != pairwiseDist.end());
            idxMin = maxMinIter.first - pairwiseDist.begin();
            //std::cout << "maxMinDist = " << maxMinDist << ", idxMin = " << idxMin << ", d = " << pairwiseDist[idxMin].first << std::endl;
        }
        size_t idxMid { (idxMin + idxMax) / 2 };
        while(idxMax > idxMin) {
            iterNum++;
            if (oracle.isMatchLess(pairwiseDist[idxMid].first)) {
                idxMax = idxMid;
            } else {
                if (idxMax - idxMin == 1)
                    idxMin++;
                else
                    idxMin = idxMid;
            }
            idxMid = (idxMin + idxMax) / 2;
        }
        return pairwiseDist[idxMid].first;
    }
    // stats
    /*
    // count number of edges
    // pairwiseDist is sorted, add edges of the same length
    int edgeNumber {idxMid};
    while(pairwiseDist[edgeNumber + 1] == pairwiseDist[edgeNumber])
        edgeNumber++;
    // add edges between diagonal points
    edgeNumber += N / 3;
    // output stats
    std::cout << idxMid << "\t" << N;
    std::cout << "\t" << iterNum;
    std::cout << "\t" << A.size() + B.size();
    std::cout << "\t" << edgeNumber << "\t";
    std::cout << (Real)(edgeNumber) / (Real)(A.size() + B.size()) << std::endl;
    */
}

// wrappers
template<class Real>
bool readDiagramPointSet(const std::string& fname, std::vector<std::pair<Real, Real>>& result)
{
    int decPrecision;
    return readDiagramPointSet(fname.c_str(), result, decPrecision);
}

template<class Real>
bool readDiagramPointSet(const char* fname, std::vector<std::pair<Real, Real>>& result)
{
    int decPrecision;
    return readDiagramPointSet(fname, result, decPrecision);
}

template<class Real>
bool readDiagramPointSet(const std::string& fname, std::vector<std::pair<Real, Real>>& result, int& decPrecision)
{
    return readDiagramPointSet(fname.c_str(), result, decPrecision);
}

// reading function
template<class Real>
bool readDiagramPointSet(const char* fname, std::vector<std::pair<Real, Real>>& result, int& decPrecision)
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
        bool fracPart = false;
        int currDecPrecision = 0;
        for(auto c : line) {
            if (c == '.') {
                fracPart = true;
            } else if (fracPart) {
                if (isdigit(c)) {
                    currDecPrecision++;
                } else {
                    fracPart = false;
                    if (currDecPrecision > decPrecision)
                        decPrecision = currDecPrecision;
                    currDecPrecision = 0;
                }
            }
        }
        Real x, y;
        std::istringstream iss(line);
        if (not(iss >> x >> y)) {
#ifndef FOR_R_TDA
            std::cerr << "Error in file " << fname << ", line number " << lineNumber << ": cannot parse \"" << line << "\"" << std::endl;
#endif
            return false;
        }
        if ( x != y ) {
            result.push_back(std::make_pair(x,y));
        } else {
#ifndef FOR_R_TDA
#ifndef VERBOSE_BOTTLENECK
            std::cerr << "Warning: in file " << fname << ", line number " << lineNumber << ", zero persistence point ignored: \"" << line << "\"" << std::endl;
#endif
#endif
        }
    }
    f.close();
    return true;
}

} // end namespace bt
} // end namespace hera
#endif // HERA_BOTTLENECK_HPP
