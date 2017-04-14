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


#include <assert.h>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <iterator>
#include <chrono>

#include "def_debug_ws.h"
#include "auction_runner_gs.h"
#include "wasserstein.h"

#ifdef FOR_R_TDA
#include "Rcpp.h"
#endif

//#define PRINT_DETAILED_TIMING

namespace geom_ws {

// *****************************
// AuctionRunnerGS 
// *****************************

AuctionRunnerGS::AuctionRunnerGS(const std::vector<DiagramPoint>& A, const std::vector<DiagramPoint>& B, const double q, const double _delta, const double _internal_p, const double _initialEpsilon, const double _epsFactor) :
    bidders(A),
    items(B),
    numBidders(A.size()),
    numItems(A.size()),
    itemsToBidders(A.size(), -1),
    biddersToItems(A.size(), -1),
    wassersteinPower(q),
    delta(_delta),
    internal_p(_internal_p),
    initialEpsilon(_initialEpsilon),
    epsilonCommonRatio(_epsFactor == 0.0 ? 5.0 : _epsFactor) 
{
    assert(initialEpsilon >= 0.0 );
    assert(epsilonCommonRatio >= 0.0 );
    assert(A.size() == B.size());
    oracle = std::unique_ptr<AuctionOracle>(new AuctionOracle(bidders, items, wassersteinPower, internal_p));
}

void AuctionRunnerGS::assignItemToBidder(IdxType itemIdx, IdxType bidderIdx)
{
    numRounds++;
    //sanityCheck();
    // only unassigned bidders should submit bids and get items
    assert(biddersToItems[bidderIdx] == -1);
    IdxType oldItemOwner = itemsToBidders[itemIdx];

    // set new owner
    biddersToItems[bidderIdx] = itemIdx;
    itemsToBidders[itemIdx] = bidderIdx;
    // remove bidder from the list of unassigned bidders
#ifdef KEEP_UNASSIGNED_ORDERED
    unassignedBidders.erase(std::make_pair(bidderIdx, bidders[bidderIdx]));
#else
    unassignedBidders.erase(bidderIdx);
#endif

    // old owner becomes unassigned
    if (oldItemOwner != -1) {
        biddersToItems[oldItemOwner] = -1;
#ifdef KEEP_UNASSIGNED_ORDERED
        unassignedBidders.insert(std::make_pair(oldItemOwner, bidders[oldItemOwner]));
#else
        unassignedBidders.insert(oldItemOwner);
#endif
    }
}


void AuctionRunnerGS::flushAssignment(void)
{
    for(auto& b2i : biddersToItems) {
        b2i = -1;
    }
    for(auto& i2b : itemsToBidders) {
        i2b = -1;
    }
    // we must flush assignment only after we got perfect matching
    assert(unassignedBidders.empty());
    // all bidders become unassigned
    for(size_t bidderIdx = 0; bidderIdx < numBidders; ++bidderIdx) {
#ifdef KEEP_UNASSIGNED_ORDERED
        unassignedBidders.insert(std::make_pair(bidderIdx, bidders[bidderIdx]));
#else
        unassignedBidders.insert(bidderIdx);
#endif
    }
    assert(unassignedBidders.size() == bidders.size());
    //oracle->adjustPrices();
}

void AuctionRunnerGS::runAuction(void)
{
    relativeError = std::numeric_limits<double>::max();
#ifdef PRINT_DETAILED_TIMING
    std::chrono::high_resolution_clock hrClock;
    std::chrono::time_point<std::chrono::high_resolution_clock> startMoment;
    startMoment = hrClock.now();
    std::vector<double> iterResults;
    std::vector<double> iterEstRelErrors;
    std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> iterTimes;
#endif
    // choose some initial epsilon
    if (initialEpsilon == 0.0)
        oracle->setEpsilon(oracle->maxVal / 4.0);
    else 
        oracle->setEpsilon(initialEpsilon);
    assert( oracle->getEpsilon() > 0 );
    int iterNum { 0 };
    bool notDone { false };
    double currentResult;
    do {
        flushAssignment();
        runAuctionPhase();
        iterNum++;
        //std::cout << "Iteration " << iterNum << " completed. " << std::endl; 
        // result is d^q
        currentResult = getDistanceToQthPowerInternal();
        double denominator = currentResult - numBidders * oracle->getEpsilon();
        currentResult = pow(currentResult, 1.0 / wassersteinPower);
#ifdef PRINT_DETAILED_TIMING
#ifndef FOR_R_TDA
        iterResults.push_back(currentResult);
        iterTimes.push_back(hrClock.now());
        std::cout << "Iteration " << iterNum << " finished. ";
        std::cout << "Current result is " << currentResult  << ", epsilon = " << oracle->getEpsilon() << std::endl;
        std::cout << "Number of rounds (cumulative): " << numRounds << std::endl;
#endif
#endif
        if ( denominator <= 0 ) {
            //std::cout << "Epsilon is too big." << std::endl;
            notDone = true;
        } else {
            denominator = pow(denominator, 1.0 / wassersteinPower);
            double numerator = currentResult - denominator;
#ifdef PRINT_DETAILED_TIMING
#ifndef FOR_R_TDA
            std::cout << " numerator: " << numerator << " denominator: " << denominator;
            std::cout << "; error bound: " << numerator / denominator << std::endl;
#endif
#endif
            relativeError = numerator / denominator;
            // if relative error is greater than delta, continue
            notDone = ( numerator / denominator > delta );
        }
        // decrease epsilon for the next iteration
        oracle->setEpsilon( oracle->getEpsilon() / epsilonCommonRatio );
        if (iterNum > maxIterNum) {
#ifndef FOR_R_TDA
            std::cerr << "Maximum iteration number exceeded, exiting. Current result is:"; 
            std::cerr << wassersteinDistance << std::endl;
#endif
            throw std::runtime_error("Maximum iteration number exceeded");
        }
    } while ( notDone );
    //printMatching();
#ifdef PRINT_DETAILED_TIMING
#ifndef FOR_R_TDA
    for(size_t iterIdx = 0; iterIdx < iterResults.size(); ++iterIdx) {
        double trueRelError = ( iterResults.at(iterIdx) - currentResult ) / currentResult;
        auto iterCumulativeTime = iterTimes.at(iterIdx) - startMoment;
        std::chrono::duration<double, std::milli> iterTime =  ( iterIdx > 0) ? iterTimes[iterIdx] - iterTimes[iterIdx - 1] : iterTimes[iterIdx] - startMoment; 
        std::cout << "iteration " << iterIdx << ", true rel. error " <<
                    trueRelError << ", elapsed time " << 
                    std::chrono::duration<double, std::milli>(iterCumulativeTime).count() << 
                    ", iteration time " << iterTime.count() << std::endl;
    }
#endif
#endif
}

void AuctionRunnerGS::runAuctionPhase(void)
{
    //std::cout << "Entered runAuctionPhase" << std::endl;
    do {
#ifdef KEEP_UNASSIGNED_ORDERED
        size_t bidderIdx = unassignedBidders.begin()->first;
#else
        size_t bidderIdx = *unassignedBidders.begin();
#endif
        auto optimalBid = oracle->getOptimalBid(bidderIdx);
        auto optimalItemIdx = optimalBid.first;
        auto bidValue = optimalBid.second;
        assignItemToBidder(optimalBid.first, bidderIdx);
        oracle->setPrice(optimalItemIdx, bidValue);
        //printDebug();
#ifdef FOR_R_TDA
        if ( numRounds % 10000 == 0 ) {
            Rcpp::checkUserInterrupt();
        }
#endif
    } while (not unassignedBidders.empty());
    //std::cout << "runAuctionPhase finished" << std::endl;

#ifdef DEBUG_AUCTION
    for(size_t bidderIdx = 0; bidderIdx < numBidders; ++bidderIdx) {
        if ( biddersToItems[bidderIdx] < 0) {
#ifndef FOR_R_TDA
            std::cerr << "After auction terminated bidder " << bidderIdx;
            std::cerr << " has no items assigned" << std::endl;
#endif
            throw std::runtime_error("Auction did not give a perfect matching");
        }
    }
#endif

}
 
double AuctionRunnerGS::getDistanceToQthPowerInternal(void)
{
    sanityCheck();
    double result = 0.0;
    //std::cout << "-------------------------------------------------------------------------\n";
    for(size_t bIdx = 0; bIdx < numBidders; ++bIdx) {
        auto pA = bidders[bIdx];
        assert( 0 <= biddersToItems[bIdx] and biddersToItems[bIdx] < static_cast<int>(items.size()) );
        auto pB = items[biddersToItems[bIdx]];
        //std::cout << "pA = " << pA << ", pB = " << pB << ", pow(distLp(pA, pB, internal_p), wassersteinPower) = " << pow(distLp(pA, pB, internal_p), wassersteinPower) << ", dist = " << distLp(pA, pB, internal_p) << std::endl;
        result += pow(distLp(pA, pB, internal_p), wassersteinPower);
    }
    //std::cout << "-------------------------------------------------------------------------\n";
    wassersteinCost = result;
    wassersteinDistance = pow(result, 1.0 / wassersteinPower);
    return result;
}

double AuctionRunnerGS::getWassersteinDistance(void)
{
    runAuction();
    return wassersteinDistance;
}

double AuctionRunnerGS::getWassersteinCost(void)
{
    runAuction();
    return wassersteinCost;
}



// Debug routines

void AuctionRunnerGS::printDebug(void)
{
#ifdef DEBUG_AUCTION
#ifndef FOR_R_TDA
    sanityCheck();
    std::cout << "**********************" << std::endl;
    std::cout << "Current assignment:" << std::endl;
    for(size_t idx = 0; idx < biddersToItems.size(); ++idx) {
        std::cout << idx << " <--> " << biddersToItems[idx] << std::endl;
    }
    std::cout << "Weights: " << std::endl;
    //for(size_t i = 0; i < numBidders; ++i) {
        //for(size_t j = 0; j < numItems; ++j) {
            //std::cout << oracle->weightMatrix[i][j] << " ";
        //}
        //std::cout << std::endl;
    //}
    std::cout << "Prices: " << std::endl;
    for(const auto price : oracle->getPrices()) {
        std::cout << price << std::endl;
    }
    std::cout << "**********************" << std::endl;
#endif
#endif
}


void AuctionRunnerGS::sanityCheck(void)
{
#ifdef DEBUG_AUCTION
    if (biddersToItems.size() != numBidders) {
#ifndef FOR_R_TDA
        std::cerr << "Wrong size of biddersToItems, must be " << numBidders << ", is " << biddersToItems.size() << std::endl;
#endif
        throw std::runtime_error("Wrong size of biddersToItems");
    }

    if (itemsToBidders.size() != numBidders) {
#ifndef FOR_R_TDA
        std::cerr << "Wrong size of itemsToBidders, must be " << numBidders << ", is " << itemsToBidders.size() << std::endl;
#endif
        throw std::runtime_error("Wrong size of itemsToBidders");
    }

    for(size_t bidderIdx = 0; bidderIdx < numBidders; ++bidderIdx) {
        if ( biddersToItems[bidderIdx] >= 0) {

            if ( std::count(biddersToItems.begin(),
                        biddersToItems.end(),
                        biddersToItems[bidderIdx]) > 1 ) {
#ifndef FOR_R_TDA
                std::cerr << "Item " << biddersToItems[bidderIdx];
                std::cerr << " appears in biddersToItems more than once" << std::endl;
#endif
                throw std::runtime_error("Duplicate in biddersToItems");
            }

            if (itemsToBidders.at(biddersToItems[bidderIdx]) != static_cast<int>(bidderIdx)) {
#ifndef FOR_R_TDA
                std::cerr << "Inconsitency: bidderIdx = " << bidderIdx;
                std::cerr << ", itemIdx in biddersToItems = ";
                std::cerr << biddersToItems[bidderIdx];
                std::cerr << ", bidderIdx in itemsToBidders = ";
                std::cerr << itemsToBidders[biddersToItems[bidderIdx]] << std::endl;
#endif
                throw std::runtime_error("inconsistent mapping");
            }
        }
    }

    for(IdxType itemIdx = 0; itemIdx < static_cast<IdxType>(numBidders); ++itemIdx) {
        if ( itemsToBidders[itemIdx] >= 0) {

            // check for uniqueness
            if ( std::count(itemsToBidders.begin(),
                        itemsToBidders.end(),
                        itemsToBidders[itemIdx]) > 1 ) {
#ifndef FOR_R_TDA
                std::cerr << "Bidder " << itemsToBidders[itemIdx];
                std::cerr << " appears in itemsToBidders more than once" << std::endl;
#endif
                throw std::runtime_error("Duplicate in itemsToBidders");
            }
            // check for consistency
            if (biddersToItems.at(itemsToBidders[itemIdx]) != static_cast<int>(itemIdx)) {
#ifndef FOR_R_TDA
                std::cerr << "Inconsitency: itemIdx = " << itemIdx;
                std::cerr << ", bidderIdx in itemsToBidders = ";
                std::cerr << itemsToBidders[itemIdx];
                std::cerr << ", itemIdx in biddersToItems= ";
                std::cerr << biddersToItems[itemsToBidders[itemIdx]] << std::endl;
#endif
                throw std::runtime_error("inconsistent mapping");
            }
        }
    }
#endif
}

void AuctionRunnerGS::printMatching(void)
{
//#ifdef DEBUG_AUCTION
#ifndef FOR_R_TDA
    sanityCheck();
    for(size_t bIdx = 0; bIdx < biddersToItems.size(); ++bIdx) {
        if (biddersToItems[bIdx] >= 0) {
            auto pA = bidders[bIdx];
            auto pB = items[biddersToItems[bIdx]];
            std::cout <<  pA << " <-> " << pB << "+" << pow(distLp(pA, pB, internal_p), wassersteinPower) << std::endl;
        } else {
            assert(false);
        }
    }
#endif
//#endif
}

} // end of namespace geom_ws
