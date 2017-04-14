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
#include <algorithm>
#include <functional>
#include <iterator>

#include "def_debug_ws.h"
#include "auction_runner_jac.h"
#include "wasserstein.h"

namespace geom_ws {

// *****************************
// AuctionRunnerJac 
// *****************************

AuctionRunnerJac::AuctionRunnerJac(const std::vector<DiagramPoint>& A, const std::vector<DiagramPoint>& B, const double q, const double _delta, const double _internal_p) :
    bidders(A),
    items(B),
    numBidders(A.size()),
    numItems(A.size()),
    itemsToBidders(A.size(), -1),
    biddersToItems(A.size(), -1),
    wassersteinPower(q),
    delta(_delta),
    internal_p(_internal_p),
    bidTable(A.size(), std::make_pair(-1, std::numeric_limits<double>::lowest()) ),
    itemReceivedBidVec(B.size(), 0 )
{
    assert(A.size() == B.size());
    oracle = std::unique_ptr<AuctionOracle>(new AuctionOracle(bidders, items, wassersteinPower, internal_p));
}

void AuctionRunnerJac::assignGoodToBidder(IdxType itemIdx, IdxType bidderIdx)
{
    //sanityCheck();
    IdxType myOldItem = biddersToItems[bidderIdx];
    IdxType currItemOwner = itemsToBidders[itemIdx];

    // set new owner
    biddersToItems[bidderIdx] = itemIdx;
    itemsToBidders[itemIdx] = bidderIdx;


    // remove bidder from the list of unassigned bidders
    unassignedBidders.erase( unassignedBiddersIterators[bidderIdx] );
    assert( 0 <= bidderIdx and bidderIdx < unassignedBiddersIterators.size() );
    unassignedBiddersIterators[bidderIdx] = unassignedBidders.end();

    if (-1 == currItemOwner) {
        // the item we want to assign does not belong to anybody,
        // just free myOldItem, if necessary
        // RE: this cannot be necessary. I submitted the best bid, hence I was
        // an unassigned bidder.
        if (myOldItem != -1) {
#ifndef FOR_R_TDA
            std::cout << "This is not happening" << std::endl;
#endif
            assert(false);
            itemsToBidders[myOldItem] = -1;
        }
    } else {
        // the current owner of itemIdx gets my old item (OK if it's -1)
        biddersToItems[currItemOwner] = myOldItem;
        // add the old owner of bids to the list of 
        if ( -1 != myOldItem ) {
#ifndef FOR_R_TDA
            std::cout << "This is not happening" << std::endl;
#endif
            assert(false);
            // if I had something, update itemsToBidders, too
            // RE: nonsense: if I had something, I am not unassigned and did not
            // sumbit any bid
            itemsToBidders[myOldItem] = currItemOwner;
        }
        unassignedBidders.push_back(currItemOwner);
        assert( unassignedBiddersIterators[currItemOwner] == unassignedBidders.end() );
        unassignedBiddersIterators[currItemOwner] = std::prev( unassignedBidders.end() );
    }
    //sanityCheck();
}


void AuctionRunnerJac::assignToBestBidder(IdxType itemIdx)
{
    assert( itemIdx >= 0 and itemIdx < static_cast<IdxType>(numItems) );
    assert( bidTable[itemIdx].first != -1);
    IdxValPair bestBid { bidTable[itemIdx] };
    assignGoodToBidder(itemIdx, bestBid.first);
    oracle->setPrice(itemIdx,  bestBid.second);
    //dynamic_cast<AuctionOracleKDTree*>(oracle)->setNai
}

void AuctionRunnerJac::clearBidTable(void)
{
    for(auto& itemWithBidIdx : itemsWithBids) {
        itemReceivedBidVec[itemWithBidIdx] = 0;
        bidTable[itemWithBidIdx].first = -1;
        bidTable[itemWithBidIdx].second = std::numeric_limits<double>::lowest();
    }
    itemsWithBids.clear();
}

void AuctionRunnerJac::submitBid(IdxType bidderIdx, const IdxValPair& itemsBidValuePair)
{
    IdxType itemIdx = itemsBidValuePair.first;
    double bidValue = itemsBidValuePair.second;
    assert( itemIdx >= 0 );
    if ( bidTable[itemIdx].second < itemsBidValuePair.second ) {
        bidTable[itemIdx].first = bidderIdx;
        bidTable[itemIdx].second = bidValue;
    }
    if (0 == itemReceivedBidVec[itemIdx]) {
        itemReceivedBidVec[itemIdx] = 1;
        itemsWithBids.push_back(itemIdx);
    }
}

void AuctionRunnerJac::printDebug(void)
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
    //std::cout << "Value matrix: " << std::endl;
    //for(size_t i = 0; i < numBidders; ++i) {
        //for(size_t j = 0; j < numItems; ++j) {
            //std::cout << oracle->weightMatrix[i][j] - oracle->prices[j] << " ";
        //}
        //std::cout << std::endl;
    //}
    std::cout << "**********************" << std::endl;
#endif
#endif
}

void AuctionRunnerJac::flushAssignment(void)
{
    for(auto& b2g : biddersToItems) {
        b2g = -1;
    }
    for(auto& g2b : itemsToBidders) {
        g2b = -1;
    }
    //oracle->adjustPrices();
}

void AuctionRunnerJac::runAuction(void)
{
    relativeError = std::numeric_limits<double>::max();
    // choose some initial epsilon
    oracle->setEpsilon(oracle->maxVal / 4.0);
    assert( oracle->getEpsilon() > 0 );
    int iterNum { 0 };
    bool notDone { false };
    do {
        flushAssignment();
        runAuctionPhase();
        iterNum++;
        //std::cout << "Iteration " << iterNum << " completed. " << std::endl; 
        // result is d^q
        double currentResult = getDistanceToQthPowerInternal();
        double denominator = currentResult - numBidders * oracle->getEpsilon();
        currentResult = pow(currentResult, 1.0 / wassersteinPower);
        //std::cout << "Current result is " << currentResult << std::endl;
        if ( denominator <= 0 ) {
            //std::cout << "Epsilon is too big." << std::endl;
            notDone = true;
        } else {
            denominator = pow(denominator, 1.0 / wassersteinPower);
            double numerator = currentResult - denominator;
            relativeError = numerator / denominator;
            //std::cout << " numerator: " << numerator << " denominator: " << denominator << std::endl;
            //std::cout << " error bound: " << numerator / denominator << std::endl;
            // if relative error is greater than delta, continue
            notDone = ( numerator / denominator > delta );
        }
        // decrease epsilon for the next iteration
        oracle->setEpsilon( oracle->getEpsilon() / epsilonCommonRatio );
        if (iterNum > maxIterNum) {
#ifndef FOR_R_TDA
            std::cerr << "Maximum iteration number exceeded, exiting. Current result is:"; 
            std::cerr << wassersteinDistance << std::endl;
            std::exit(1);
#else
            throw std::runtime_error("Max. iter. exceeded");
#endif
        }
    } while ( notDone );
    //printMatching();
}

void AuctionRunnerJac::runAuctionPhase(void)
{
    //std::cout << "Entered runAuctionPhase" << std::endl;
    //int numUnassignedBidders { 0 };

    // at the beginning of a phase all bidders are unassigned
    unassignedBidders.clear();
    unassignedBiddersIterators.clear();
    for(size_t bidderIdx = 0; bidderIdx < numBidders; ++bidderIdx) {
        unassignedBidders.push_back(bidderIdx);
        unassignedBiddersIterators.push_back( std::prev( unassignedBidders.end() ));
    }
    do {
        // bidding phase
        clearBidTable();
        for(const auto bidderIdx : unassignedBidders) {
            submitBid(bidderIdx, oracle->getOptimalBid(bidderIdx));
        }
        //std::cout << "Number of unassignedBidders: " << unassignedBidders.size() << std::endl;

        // assignment phase
        // todo: maintain list of items that received a bid
        for(auto itemIdx : itemsWithBids ) {
            assignToBestBidder(itemIdx);
        }
        //std::cout << "Assignment phase done" << std::endl;
        //sanityCheck();
        //printDebug();
    } while (unassignedBidders.size() > 0);
    //std::cout << "runAuctionPhase finished" << std::endl;


#ifdef DEBUG_AUCTION
    for(size_t bidderIdx = 0; bidderIdx < numBidders; ++bidderIdx) {
        if ( biddersToItems[bidderIdx] < 0) {
            std::cerr << "After auction terminated bidder " << bidderIdx;
            std::cerr << " has no items assigned" << std::endl;
            throw "Auction did not give a perfect matching";
        }
    }
#endif

}
 
// assertion: the matching must be perfect
double AuctionRunnerJac::getDistanceToQthPowerInternal(void)
{
    sanityCheck();
    double result = 0.0;
    for(size_t bIdx = 0; bIdx < numBidders; ++bIdx) {
        auto pA = bidders[bIdx];
        assert( 0 <= biddersToItems[bIdx] and biddersToItems[bIdx] < items.size() );
        auto pB = items[biddersToItems[bIdx]];
        result += pow(distLp(pA, pB, internal_p), wassersteinPower);
    }
    wassersteinCost = result;
    wassersteinDistance = pow(result, 1.0 / wassersteinPower);
    return result;
}

double AuctionRunnerJac::getWassersteinDistance(void)
{
    runAuction();
    return wassersteinDistance;
}

double AuctionRunnerJac::getWassersteinCost(void)
{
    runAuction();
    return wassersteinCost;
}



void AuctionRunnerJac::sanityCheck(void)
{
#ifdef DEBUG_AUCTION
#ifndef FOR_R_TDA
    if (biddersToItems.size() != numBidders) {
        std::cerr << "Wrong size of biddersToItems, must be " << numBidders << ", is " << biddersToItems.size() << std::endl;
        throw "Wrong size of biddersToItems";
    }

    if (itemsToBidders.size() != numBidders) {
        std::cerr << "Wrong size of itemsToBidders, must be " << numBidders << ", is " << itemsToBidders.size() << std::endl;
        throw "Wrong size of itemsToBidders";
    }

    for(size_t bidderIdx = 0; bidderIdx < numBidders; ++bidderIdx) {
        if ( biddersToItems[bidderIdx] >= 0) {

            if ( std::count(biddersToItems.begin(),
                        biddersToItems.end(),
                        biddersToItems[bidderIdx]) > 1 ) {
                std::cerr << "Good " << biddersToItems[bidderIdx];
                std::cerr << " appears in biddersToItems more than once" << std::endl;
                throw "Duplicate in biddersToItems";
            }

            if (itemsToBidders.at(biddersToItems[bidderIdx]) != static_cast<int>(bidderIdx)) {
                std::cerr << "Inconsitency: bidderIdx = " << bidderIdx;
                std::cerr << ", itemIdx in biddersToItems = ";
                std::cerr << biddersToItems[bidderIdx];
                std::cerr << ", bidderIdx in itemsToBidders = ";
                std::cerr << itemsToBidders[biddersToItems[bidderIdx]] << std::endl;
                throw "inconsistent mapping";
            }
        }
    }

    for(IdxType itemIdx = 0; itemIdx < static_cast<IdxType>(numBidders); ++itemIdx) {
        if ( itemsToBidders[itemIdx] >= 0) {

            // check for uniqueness
            if ( std::count(itemsToBidders.begin(),
                        itemsToBidders.end(),
                        itemsToBidders[itemIdx]) > 1 ) {
                std::cerr << "Bidder " << itemsToBidders[itemIdx];
                std::cerr << " appears in itemsToBidders more than once" << std::endl;
                throw "Duplicate in itemsToBidders";
            }
            // check for consistency
            if (biddersToItems.at(itemsToBidders[itemIdx]) != static_cast<int>(itemIdx)) {
                std::cerr << "Inconsitency: itemIdx = " << itemIdx;
                std::cerr << ", bidderIdx in itemsToBidders = ";
                std::cerr << itemsToBidders[itemIdx];
                std::cerr << ", itemIdx in biddersToItems= ";
                std::cerr << biddersToItems[itemsToBidders[itemIdx]] << std::endl;
                throw "inconsistent mapping";
            }
        }
    }
#endif
#endif
}

void AuctionRunnerJac::printMatching(void)
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
