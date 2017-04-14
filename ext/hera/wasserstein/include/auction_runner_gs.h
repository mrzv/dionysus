/*
 
Copyright (c) 2016, M. Kerber, D. Morozov, A. Nigmetov
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

#ifndef AUCTION_RUNNER_GS_H
#define AUCTION_RUNNER_GS_H

#include <unordered_set>

#include "auction_oracle.h"

//#define KEEP_UNASSIGNED_ORDERED
// if this symbol is defined, 
// unassigned bidders are processed in a lexicographic order.
// See LexicogrCompDiagramPoint comparator.


namespace geom_ws {

//using AuctionOracle = AuctionOracleLazyHeapRestricted;
using AuctionOracle = AuctionOracleKDTreeRestricted;

#ifdef KEEP_UNASSIGNED_ORDERED
using IdxPointPair = std::pair<size_t, DiagramPoint>;

struct LexicogrCompDiagramPoint {
    bool operator ()(const IdxPointPair& a, const IdxPointPair& b) {
        const auto& p1 = a.second;
        const auto& p2 = b.second;

        return ( (not p1.isDiagonal() and p2.isDiagonal()) or
                ( p1.isDiagonal() == p2.isDiagonal() and p1.getRealX() < p2.getRealX() ) or
                ( p1.isDiagonal() == p2.isDiagonal() and p1.getRealX() == p2.getRealX() and p1.getRealY() < p2.getRealY() ) or
                ( p1.isDiagonal() == p2.isDiagonal() and p1.getRealX() == p2.getRealX() and p1.getRealY() == p2.getRealY() and a.first < b.first ) );
    }
};

using OrderedUnassignedKeeper = std::set<IdxPointPair, LexicogrCompDiagramPoint>;
#endif

// the two parameters that you can tweak in auction algorithm are:
// 1. epsilonCommonRatio
// 2. maxIterNum

class AuctionRunnerGS {
public:
    AuctionRunnerGS(const std::vector<DiagramPoint>& A,
                    const std::vector<DiagramPoint>& B,
                    const double q,
                    const double _delta,
                    const double _internal_p,
                    const double _initialEpsilon,
                    const double _epsFactor);
    void setEpsilon(double newVal) { assert(epsilon > 0.0); epsilon = newVal; };
    double getEpsilon() const { return epsilon; }
    double getWassersteinDistance();
    double getWassersteinCost();
    double getRelativeError() const { return relativeError; };
    static constexpr int maxIterNum { 25 }; // maximal number of iterations of epsilon-scaling
private:
    // private data
    std::vector<DiagramPoint> bidders, items;
    const size_t numBidders;
    const size_t numItems;
    std::vector<IdxType> itemsToBidders;
    std::vector<IdxType> biddersToItems;
    double wassersteinPower;
    double epsilon;
    double delta;
    double internal_p;
    double initialEpsilon;
    double epsilonCommonRatio; // next epsilon = current epsilon / epsilonCommonRatio
    double weightAdjConst;
    double wassersteinDistance;
    double wassersteinCost;
    double relativeError;
    // to get the 2 best items
    std::unique_ptr<AuctionOracle> oracle;
#ifdef KEEP_UNASSIGNED_ORDERED
    OrderedUnassignedKeeper unassignedBidders;
#else
    std::unordered_set<size_t> unassignedBidders;
#endif
    // private methods
    void assignItemToBidder(const IdxType bidderIdx, const IdxType itemsIdx);
    void clearBidTable(void);
    void runAuction(void);
    void runAuctionPhase(void);
    void flushAssignment(void);

    // for debug only
    void sanityCheck(void);
    void printDebug(void);
    int countUnhappy(void);
    void printMatching(void);
    double getDistanceToQthPowerInternal(void);
    int numRounds { 0 };
};

} // end of namespace geom_ws

#endif
