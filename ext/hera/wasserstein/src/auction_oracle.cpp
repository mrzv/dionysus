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
#include "auction_oracle.h"

namespace geom_ws {

AuctionOracleAbstract::AuctionOracleAbstract(const std::vector<DiagramPoint>& _bidders, const std::vector<DiagramPoint>& _items, const double _wassersteinPower, const double _internal_p) :
    bidders(_bidders),
    items(_items),
    prices(items.size(), 0.0),
    wassersteinPower(_wassersteinPower),
    internal_p(_internal_p)
{
}

double AuctionOracleAbstract::getValueForBidder(size_t bidderIdx, size_t itemIdx)
{
    return pow(distLp(bidders[bidderIdx], items[itemIdx], internal_p), wassersteinPower) + prices[itemIdx];
}

// *****************************
// AuctionOracleLazyHeap
// *****************************

AuctionOracleLazyHeap::AuctionOracleLazyHeap(const std::vector<DiagramPoint>& b, 
                                     const std::vector<DiagramPoint>& g, 
                                     double _wassersteinPower,
                                     const double internal_p) :
    AuctionOracleAbstract(b, g, _wassersteinPower, internal_p),
    maxVal(std::numeric_limits<double>::min()),
    biddersUpdateMoments(b.size(), 0),
    updateCounter(0)
{
    assert(b.size() == g.size() );
    assert(b.size() > 1);

    weightMatrix.reserve(b.size());
    //const double maxDistUpperBound = 3 * getFurthestDistance3Approx(b, g);
    //weightAdjConst = pow(maxDistUpperBound, wassersteinPower);
    // init weight matrix
    for(const auto& pointA : bidders) {
        std::vector<double> weightVec;
        weightVec.clear();
        weightVec.reserve(b.size());
        for(const auto& pointB : items) {
            double val = pow(distLp(pointA, pointB, internal_p), wassersteinPower);
            if (val > maxVal) {
                maxVal = val;
            }
            weightVec.push_back( val );
        }
        weightMatrix.push_back(weightVec);
    }
    fillInLossesHeap();
    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        updateList.push_back(std::make_pair(static_cast<IdxType>(itemIdx), 0));
    }
    for(auto updateListIter = updateList.begin(); updateListIter != updateList.end(); ++updateListIter) {
        itemsIterators.push_back(updateListIter);
    }
}

void AuctionOracleLazyHeap::updateQueueForBidder(IdxType bidderIdx)
{
    assert(0 <= bidderIdx and bidderIdx < static_cast<int>(bidders.size()));
    assert(bidderIdx < static_cast<int>(biddersUpdateMoments.size()));

    int bidderLastUpdateTime = biddersUpdateMoments[bidderIdx];
    auto iter = updateList.begin();
    while (iter != updateList.end() and iter->second >= bidderLastUpdateTime) {
        IdxType itemIdx = iter->first;
        IdxValPair newVal { itemIdx, weightMatrix[bidderIdx][itemIdx] + prices[itemIdx]};
        // to-do: change indexing of lossesHeapHandles
        lossesHeap[bidderIdx]->decrease(lossesHeapHandles[bidderIdx][itemIdx], newVal);
        iter++;
    }
    biddersUpdateMoments[bidderIdx] = updateCounter;
}

void AuctionOracleLazyHeap::fillInLossesHeap(void)
{
    for(size_t bidderIdx = 0; bidderIdx < bidders.size(); ++bidderIdx) {
        lossesHeap.push_back( new LossesHeap() );
        std::vector<LossesHeap::handle_type> handlesVec;
        lossesHeapHandles.push_back(handlesVec);
        for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
            IdxValPair vp { itemIdx, weightMatrix[bidderIdx][itemIdx] + prices[itemIdx] };
            lossesHeapHandles[bidderIdx].push_back(  lossesHeap[bidderIdx]->push(vp) );
        }
    }
}

AuctionOracleLazyHeap::~AuctionOracleLazyHeap()
{
    for(auto h : lossesHeap) {
        delete h;
    }
}

void AuctionOracleLazyHeap::setPrice(IdxType itemIdx, double newPrice)
{
    assert( prices.at(itemIdx) < newPrice );
#ifdef DEBUG_AUCTION
#ifndef FOR_R_TDA
    std::cout << "price incremented by " <<  prices.at(itemIdx) - newPrice << std::endl;
#endif
#endif
    prices[itemIdx] = newPrice;
    // lazy: record the moment we updated the price of the items,
    // do not update queues.
    // 1. move the items with updated price to the front of the updateList,
    updateList.splice(updateList.begin(), updateList, itemsIterators[itemIdx]);
    // 2. record the moment we updated the price and increase the counter
    updateList.front().second = updateCounter++;
}

// subtract min. price from all prices
void AuctionOracleLazyHeap::adjustPrices(void)
{
    double minPrice = *(std::min_element(prices.begin(), prices.end()));
    std::transform(prices.begin(), prices.end(), prices.begin(), [minPrice](double a) { return a - minPrice; });
}

DebugOptimalBid AuctionOracleLazyHeap::getOptimalBidDebug(IdxType bidderIdx)
{
    assert(bidderIdx >=0 and bidderIdx < static_cast<IdxType>(bidders.size()) );
    assert(lossesHeap.at(bidderIdx) != nullptr);
    assert(lossesHeap[bidderIdx]->size() >= 2);

    updateQueueForBidder(bidderIdx);
    DebugOptimalBid result;
    
    auto pHeap = lossesHeap[bidderIdx];
    auto topIter = pHeap->ordered_begin(); 
    result.bestItemIdx = topIter->first;
    result.bestItemValue = topIter->second;
    ++topIter; // now points to the second-best items
    result.secondBestItemValue = topIter->second;
    result.secondBestItemIdx = topIter->first;

#ifndef FOR_R_TDA
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemValue = " << bestItemValue << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestValue = " << secondBestItemValue << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemsDist= " << (weightAdjConst -  bestItemValue) << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestDist= " << (weightAdjConst - secondBestItemValue) << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
#endif

    return result;
}

IdxValPair AuctionOracleLazyHeap::getOptimalBid(const IdxType bidderIdx) 
{
    assert(bidderIdx >=0 and bidderIdx < static_cast<IdxType>(bidders.size()) );
    assert(lossesHeap.at(bidderIdx) != nullptr);
    assert(lossesHeap[bidderIdx]->size() >= 2);

    updateQueueForBidder(bidderIdx);
    
    auto pHeap = lossesHeap[bidderIdx];
    auto topIter = pHeap->ordered_begin(); 
    IdxType bestItemIdx = topIter->first;
    double bestItemValue = topIter->second;
    ++topIter; // now points to the second-best items
    double secondBestItemValue = topIter->second;

#ifndef FOR_R_TDA
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemValue = " << bestItemValue << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestValue = " << secondBestItemValue << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemsDist= " << (weightAdjConst -  bestItemValue) << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestDist= " << (weightAdjConst - secondBestItemValue) << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
#endif

    // bid value: price + value difference + epsilon
    return std::make_pair(bestItemIdx, 
                          prices[bestItemIdx] + 
                          ( secondBestItemValue - bestItemValue ) +
                          epsilon );
}

// *****************************
// AuctionOracleLazyHeapRestricted
// *****************************

AuctionOracleLazyHeapRestricted::AuctionOracleLazyHeapRestricted(const std::vector<DiagramPoint>& b, 
                                     const std::vector<DiagramPoint>& g, 
                                     double _wassersteinPower,
                                     double internal_p) :
    AuctionOracleAbstract(b, g, _wassersteinPower),
    maxVal(std::numeric_limits<double>::min()),
    biddersUpdateMoments(b.size(), 0),
    updateCounter(0),
    heapHandlesIndices(items.size(), std::numeric_limits<size_t>::max()),
    bestDiagonalItemsComputed(false)
{
    assert(b.size() == g.size() );
    assert(b.size() > 1);

    weightMatrix.reserve(b.size());
    //const double maxDistUpperBound = 3 * getFurthestDistance3Approx(b, g);
    //weightAdjConst = pow(maxDistUpperBound, wassersteinPower);
    // init weight matrix
    for(const auto& pointA : bidders) {
        std::vector<double> weightVec;
        weightVec.clear();
        weightVec.reserve(b.size());
        for(const auto& pointB : items) {
            double val = pow(distLp(pointA, pointB, internal_p), wassersteinPower);
            weightVec.push_back( val );
            if ( val > maxVal )
                maxVal = val;
        }
        weightMatrix.push_back(weightVec);
    }
    fillInLossesHeap();
    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        updateList.push_back(std::make_pair(static_cast<IdxType>(itemIdx), 0));
    }
    for(auto updateListIter = updateList.begin(); updateListIter != updateList.end(); ++updateListIter) {
        itemsIterators.push_back(updateListIter);
    }

    size_t handleIdx {0};
    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        if (items[itemIdx].isDiagonal() ) {
            heapHandlesIndices[itemIdx] = handleIdx++;
            diagHeapHandles.push_back(diagItemsHeap.push(std::make_pair(itemIdx, 0)));
        }
    }
     // todo: this must be done in readFiles procedure
}

void AuctionOracleLazyHeapRestricted::updateQueueForBidder(IdxType bidderIdx)
{
    assert(0 <= bidderIdx and bidderIdx < static_cast<int>(bidders.size()));
    assert(bidderIdx < static_cast<int>(biddersUpdateMoments.size()));
    assert(lossesHeap[bidderIdx] != nullptr );

    int bidderLastUpdateTime = biddersUpdateMoments[bidderIdx];
    auto iter = updateList.begin();
    while (iter != updateList.end() and iter->second >= bidderLastUpdateTime) {
        IdxType itemIdx = iter->first;
        size_t handleIdx = itemsIndicesForHeapHandles[bidderIdx][itemIdx];
        if (handleIdx  < items.size() ) {
            IdxValPair newVal { itemIdx, weightMatrix[bidderIdx][itemIdx] + prices[itemIdx]};
            // to-do: change indexing of lossesHeapHandles
            lossesHeap[bidderIdx]->decrease(lossesHeapHandles[bidderIdx][handleIdx], newVal);
        }
        iter++;
    }
    biddersUpdateMoments[bidderIdx] = updateCounter;
}

void AuctionOracleLazyHeapRestricted::fillInLossesHeap(void)
{
    for(size_t bidderIdx = 0; bidderIdx < bidders.size(); ++bidderIdx) {
        DiagramPoint bidder { bidders[bidderIdx] };
        // no heap for diagonal bidders
        if ( bidder.isDiagonal() ) {
            lossesHeap.push_back( nullptr );
            lossesHeapHandles.push_back(std::vector<LossesHeap::handle_type>());
            itemsIndicesForHeapHandles.push_back( std::vector<size_t>() );
            continue;
        } else {
            lossesHeap.push_back( new LossesHeap() );
            assert( lossesHeap.at(bidderIdx) != nullptr );
            itemsIndicesForHeapHandles.push_back( std::vector<size_t>(items.size(), std::numeric_limits<size_t>::max() ) );

            std::vector<LossesHeap::handle_type> handlesVec;
            lossesHeapHandles.push_back(handlesVec);
            size_t handleIdx { 0 };
            for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
                assert( itemsIndicesForHeapHandles.at(bidderIdx).at(itemIdx) > 0 );
                DiagramPoint item { items[itemIdx] };
                if ( item.isNormal() ) {
                    // item can be assigned to bidder, store in heap 
                    IdxValPair vp { itemIdx, weightMatrix[bidderIdx][itemIdx] + prices[itemIdx] };
                    lossesHeapHandles[bidderIdx].push_back(  lossesHeap[bidderIdx]->push(vp) );
                    // keep corresponding index in itemsIndicesForHeapHandles
                    itemsIndicesForHeapHandles[bidderIdx][itemIdx] = handleIdx++;
                }
            }
        }
    }
}

AuctionOracleLazyHeapRestricted::~AuctionOracleLazyHeapRestricted()
{
    for(auto h : lossesHeap) {
        delete h;
    }
}

void AuctionOracleLazyHeapRestricted::setPrice(IdxType itemIdx, double newPrice)
{
    assert( prices.at(itemIdx) < newPrice );
#ifdef DEBUG_AUCTION
#ifndef FOR_R_TDA
    std::cout << "price incremented by " <<  prices.at(itemIdx) - newPrice << std::endl;
#endif
#endif
    prices[itemIdx] = newPrice;
    if (items[itemIdx].isNormal() ) {
        // lazy: record the moment we updated the price of the items,
        // do not update queues.
        // 1. move the items with updated price to the front of the updateList,
        updateList.splice(updateList.begin(), updateList, itemsIterators[itemIdx]);
        // 2. record the moment we updated the price and increase the counter
        updateList.front().second = updateCounter++;
    } else {
        // diagonal items are stored in one heap
        diagItemsHeap.decrease(diagHeapHandles[heapHandlesIndices[itemIdx]], std::make_pair(itemIdx, newPrice));
        bestDiagonalItemsComputed = false;
    }
}

// subtract min. price from all prices
void AuctionOracleLazyHeapRestricted::adjustPrices(void)
{
}

DebugOptimalBid AuctionOracleLazyHeapRestricted::getOptimalBidDebug(IdxType bidderIdx)
{
    DebugOptimalBid result;
    assert(bidderIdx >=0 and bidderIdx < static_cast<IdxType>(bidders.size()) );

    DiagramPoint bidder = bidders[bidderIdx];
    std::vector<IdxValPair> candItems;
    // corresponding point is always considered as a candidate

    size_t projItemIdx = bidderIdx;
    assert( 0 <= projItemIdx and projItemIdx < items.size() );
    DiagramPoint projItem = items[projItemIdx];
    assert(projItem.type != bidder.type);
    //assert(projItem.projId == bidder.id);
    //assert(projItem.id == bidder.projId);
    // todo: store precomputed distance?
    double projItemValue = pow(distLp(bidder, projItem, internal_p), wassersteinPower) + prices[projItemIdx];
    candItems.push_back( std::make_pair(projItemIdx, projItemValue) );
 
    if (bidder.isNormal()) {
        assert(lossesHeap.at(bidderIdx) != nullptr);
        assert(lossesHeap[bidderIdx]->size() >= 2);
        updateQueueForBidder(bidderIdx);
        auto pHeap = lossesHeap[bidderIdx];
        assert( pHeap != nullptr );
        auto topIter = pHeap->ordered_begin(); 
        candItems.push_back( *topIter );
        ++topIter; // now points to the second-best items
        candItems.push_back( *topIter );
        std::sort(candItems.begin(), candItems.end(), CompPairsBySecondStruct());
        assert(candItems[1].second >= candItems[0].second);
    } else {
        // for diagonal bidder the only normal point has already been added
        // the other 2 candidates are diagonal items only, get from the heap
        // with prices
        assert(diagItemsHeap.size() > 1);
        auto topDiagIter = diagItemsHeap.ordered_begin();
        auto topDiag1 = *topDiagIter++;
        auto topDiag2 = *topDiagIter;
        candItems.push_back(topDiag1);
        candItems.push_back(topDiag2);
        std::sort(candItems.begin(), candItems.end(), CompPairsBySecondStruct());
        assert(candItems.size() == 3);
        assert(candItems[2].second >= candItems[1].second);
        assert(candItems[1].second >= candItems[0].second);
    }
    
    result.bestItemIdx = candItems[0].first;
    result.secondBestItemIdx = candItems[1].first;
    result.bestItemValue = candItems[0].second;
    result.secondBestItemValue = candItems[1].second;

    // checking code

    //DebugOptimalBid debugMyResult(result);
    //DebugOptimalBid debugNaiveResult;
    //debugNaiveResult.bestItemValue = 1e20;
    //debugNaiveResult.secondBestItemValue = 1e20;
    //double currItemValue;
    //for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        //if ( bidders[bidderIdx].type != items[itemIdx].type and
                //bidders[bidderIdx].projId != items[itemIdx].id)
            //continue;

        //currItemValue = pow(distLp(bidders[bidderIdx], items[itemIdx]), wassersteinPower) + prices[itemIdx];
        //if (currItemValue < debugNaiveResult.bestItemValue) {
            //debugNaiveResult.bestItemValue = currItemValue;
            //debugNaiveResult.bestItemIdx  = itemIdx;
        //}
    //}

    //for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        //if (itemIdx == debugNaiveResult.bestItemIdx) {
            //continue;
        //}
        //if ( bidders[bidderIdx].type != items[itemIdx].type and
                //bidders[bidderIdx].projId != items[itemIdx].id)
            //continue;

        //currItemValue = pow(distLp(bidders[bidderIdx], items[itemIdx]), wassersteinPower) + prices[itemIdx];
        //if (currItemValue < debugNaiveResult.secondBestItemValue) {
            //debugNaiveResult.secondBestItemValue = currItemValue;
            //debugNaiveResult.secondBestItemIdx = itemIdx;
        //}
    //}

    //if ( fabs( debugMyResult.bestItemValue - debugNaiveResult.bestItemValue ) > 1e-6 or
            //fabs( debugNaiveResult.secondBestItemValue - debugMyResult.secondBestItemValue) > 1e-6 ) {
        //std::cerr << "bidderIdx = " << bidderIdx << "; ";
        //std::cerr << bidders[bidderIdx] << std::endl;
        //for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
            //std::cout << itemIdx << ": " << items[itemIdx] << "; price = " << prices[itemIdx] << std::endl;
        //}
        //std::cerr << "debugMyResult: " << debugMyResult << std::endl;
        //std::cerr << "debugNaiveResult: " << debugNaiveResult << std::endl;
        //auto pHeap = lossesHeap[bidderIdx];
        //assert( pHeap != nullptr );
        //for(auto topIter = pHeap->ordered_begin(); topIter != pHeap->ordered_end(); ++topIter) {
            //std::cerr << "in heap: " << topIter->first << ": " << topIter->second << "; real value = " << distLp(bidder, items[topIter->first]) + prices[topIter->first] << std::endl;
        //}
        //for(auto ci : candItems) {
            //std::cout << "ci.idx = " << ci.first << ", value = " << ci.second << std::endl;
        //}

        ////std::cerr << "twoBestItems: " << twoBestItems[0].d << " " << twoBestItems[1].d << std::endl;
        //assert(false);
    //}


    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemValue = " << bestItemValue << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestValue = " << secondBestItemValue << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemsDist= " << (weightAdjConst -  bestItemValue) << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestDist= " << (weightAdjConst - secondBestItemValue) << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;

    return result;
}

IdxValPair AuctionOracleLazyHeapRestricted::getOptimalBid(const IdxType bidderIdx) 
{
    IdxType bestItemIdx;
    //IdxType secondBestItemIdx;
    double bestItemValue;
    double secondBestItemValue;

    auto& bidder = bidders[bidderIdx];
    IdxType projItemIdx = bidderIdx;
    assert( 0 <= projItemIdx and projItemIdx < items.size() );
    DiagramPoint projItem = items[projItemIdx];
    assert(projItem.type != bidder.type);
    //assert(projItem.projId == bidder.id);
    //assert(projItem.id == bidder.projId);
    // todo: store precomputed distance?
    double projItemValue = pow(distLp(bidder, projItem, internal_p), wassersteinPower) + prices[projItemIdx];
   
    if (bidder.isDiagonal()) {
        // for diagonal bidder the only normal point has already been added
        // the other 2 candidates are diagonal items only, get from the heap
        // with prices
        assert(diagItemsHeap.size() > 1);
        if (!bestDiagonalItemsComputed) {
            auto topDiagIter = diagItemsHeap.ordered_begin();
            bestDiagonalItemIdx = topDiagIter->first;
            bestDiagonalItemValue = topDiagIter->second;
            topDiagIter++;
            secondBestDiagonalItemIdx = topDiagIter->first;
            secondBestDiagonalItemValue = topDiagIter->second;
            bestDiagonalItemsComputed = true;
        }

        if ( projItemValue < bestDiagonalItemValue) {
            bestItemIdx = projItemIdx;
            bestItemValue = projItemValue;
            secondBestItemValue = bestDiagonalItemValue;
            //secondBestItemIdx = bestDiagonalItemIdx;
        } else if (projItemValue < secondBestDiagonalItemValue) {
            bestItemIdx = bestDiagonalItemIdx;
            bestItemValue = bestDiagonalItemValue;
            secondBestItemValue = projItemValue;
            //secondBestItemIdx = projItemIdx;
        } else {
            bestItemIdx = bestDiagonalItemIdx;
            bestItemValue = bestDiagonalItemValue;
            secondBestItemValue = secondBestDiagonalItemValue;
            //secondBestItemIdx = secondBestDiagonalItemIdx;
        }
    } else {
        // for normal bidder get 2 best items among non-diagonal (=normal) points 
        // from the corresponding heap 
        assert(diagItemsHeap.size() > 1);
        updateQueueForBidder(bidderIdx);
        auto topNormIter = lossesHeap[bidderIdx]->ordered_begin();
        IdxType bestNormalItemIdx { topNormIter->first };
        double bestNormalItemValue { topNormIter->second };
        topNormIter++;
        double secondBestNormalItemValue { topNormIter->second };
        //IdxType secondBestNormalItemIdx { topNormIter->first };
 
        if ( projItemValue < bestNormalItemValue) {
            bestItemIdx = projItemIdx;
            bestItemValue = projItemValue;
            secondBestItemValue = bestNormalItemValue;
            //secondBestItemIdx = bestNormalItemIdx;
        } else if (projItemValue < secondBestNormalItemValue) {
            bestItemIdx = bestNormalItemIdx;
            bestItemValue = bestNormalItemValue;
            secondBestItemValue = projItemValue;
            //secondBestItemIdx = projItemIdx;
        } else {
            bestItemIdx = bestNormalItemIdx;
            bestItemValue = bestNormalItemValue;
            secondBestItemValue = secondBestNormalItemValue;
            //secondBestItemIdx = secondBestNormalItemIdx;
        }
    }

    IdxValPair result;

    assert( secondBestItemValue >= bestItemValue );

    result.first = bestItemIdx;
    result.second = ( secondBestItemValue - bestItemValue ) + prices[bestItemIdx] + epsilon;


    // checking code

    //DebugOptimalBid debugMyResult;
    //debugMyResult.bestItemIdx = bestItemIdx;
    //debugMyResult.bestItemValue = bestItemValue;
    //debugMyResult.secondBestItemIdx = secondBestItemIdx;
    //debugMyResult.secondBestItemValue = secondBestItemValue;
    //DebugOptimalBid debugNaiveResult;
    //debugNaiveResult.bestItemValue = 1e20;
    //debugNaiveResult.secondBestItemValue = 1e20;
    //double currItemValue;
    //for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        //if ( bidders[bidderIdx].type != items[itemIdx].type and
                //bidders[bidderIdx].projId != items[itemIdx].id)
            //continue;

        //currItemValue = pow(distLp(bidders[bidderIdx], items[itemIdx]), wassersteinPower) + prices[itemIdx];
        //if (currItemValue < debugNaiveResult.bestItemValue) {
            //debugNaiveResult.bestItemValue = currItemValue;
            //debugNaiveResult.bestItemIdx  = itemIdx;
        //}
    //}

    //for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        //if (itemIdx == debugNaiveResult.bestItemIdx) {
            //continue;
        //}
        //if ( bidders[bidderIdx].type != items[itemIdx].type and
                //bidders[bidderIdx].projId != items[itemIdx].id)
            //continue;

        //currItemValue = pow(distLp(bidders[bidderIdx], items[itemIdx]), wassersteinPower) + prices[itemIdx];
        //if (currItemValue < debugNaiveResult.secondBestItemValue) {
            //debugNaiveResult.secondBestItemValue = currItemValue;
            //debugNaiveResult.secondBestItemIdx = itemIdx;
        //}
    //}
    ////std::cout << "got naive result" << std::endl;

    //if ( fabs( debugMyResult.bestItemValue - debugNaiveResult.bestItemValue ) > 1e-6 or
            //fabs( debugNaiveResult.secondBestItemValue - debugMyResult.secondBestItemValue) > 1e-6 ) {
        //std::cerr << "bidderIdx = " << bidderIdx << "; ";
        //std::cerr << bidders[bidderIdx] << std::endl;
        //for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
            //std::cout << itemIdx << ": " << items[itemIdx] << "; price = " << prices[itemIdx] << std::endl;
        //}
        //std::cerr << "debugMyResult: " << debugMyResult << std::endl;
        //std::cerr << "debugNaiveResult: " << debugNaiveResult << std::endl;
        //auto pHeap = lossesHeap[bidderIdx];
        //if ( pHeap != nullptr ) {
            //for(auto topIter = pHeap->ordered_begin(); topIter != pHeap->ordered_end(); ++topIter) {
                //std::cerr << "in heap: " << topIter->first << ": " << topIter->second << "; real value = " << distLp(bidder, items[topIter->first]) + prices[topIter->first] << std::endl;
            //}
        //}
        ////for(auto ci : candItems) {
            ////std::cout << "ci.idx = " << ci.first << ", value = " << ci.second << std::endl;
        ////}

        ////std::cerr << "twoBestItems: " << twoBestItems[0].d << " " << twoBestItems[1].d << std::endl;
        //assert(false);
    // }
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemValue = " << bestItemValue << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestValue = " << secondBestItemValue << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemsDist= " << (weightAdjConst -  bestItemValue) << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestDist= " << (weightAdjConst - secondBestItemValue) << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;

    return result;
}


// *****************************
// AuctionOracleKDTree
// *****************************

AuctionOracleKDTree::AuctionOracleKDTree(const std::vector<DiagramPoint>& _bidders, 
        const std::vector<DiagramPoint>& _items, 
        double _wassersteinPower,
        double internal_p) :
    AuctionOracleAbstract(_bidders, _items, _wassersteinPower, internal_p),
    heapHandlesIndices(items.size(), std::numeric_limits<size_t>::max()),
    kdtreeItems(items.size(), std::numeric_limits<size_t>::max())
{
    //assert(wassersteinPower == 1.0); // temporarily, to-do: update dnn to search with any q
    size_t dnnItemIdx { 0 };
    size_t trueIdx { 0 };
    dnnPoints.clear();
    // store normal items in kd-tree
    for(const auto& g : items) {
        if (g.isNormal()) {
            kdtreeItems[trueIdx] = dnnItemIdx;
            // index of items is id of dnn-point
            DnnPoint p(trueIdx);
            p[0] = g.getRealX();
            p[1] = g.getRealY();
            dnnPoints.push_back(p);
            assert(dnnItemIdx == dnnPoints.size() - 1);
            dnnItemIdx++;
        }
        trueIdx++;
    }

    assert(dnnPoints.size() < items.size() );
    for(size_t i = 0; i < dnnPoints.size(); ++i) {
        dnnPointHandles.push_back(&dnnPoints[i]);
    }
    DnnTraits traits;
    traits.internal_p = internal_p;
    kdtree = new dnn::KDTree<DnnTraits>(traits, dnnPointHandles, wassersteinPower);

    size_t dnnItemIdxAll { 0 };
    dnnPointsAll.clear();
    // store all items in kd-tree
    for(const auto& g : items) {
        DnnPoint p(dnnItemIdxAll++);
        p[0] = g.getRealX();
        p[1] = g.getRealY();
        dnnPointsAll.push_back(p);
        assert(dnnItemIdxAll == dnnPointsAll.size());
    }

    for(size_t i = 0; i < dnnPointsAll.size(); ++i) {
        dnnPointHandlesAll.push_back(&dnnPointsAll[i]);
    }
    kdtreeAll = new dnn::KDTree<DnnTraits>(traits, dnnPointHandlesAll, wassersteinPower);
    
    size_t handleIdx {0};
    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        if (items[itemIdx].isDiagonal() ) {
            heapHandlesIndices[itemIdx] = handleIdx++;
            diagHeapHandles.push_back(diagItemsHeap.push(std::make_pair(itemIdx, 0)));
        }
    }
    //to-do: remove maxVal from 
    maxVal = 3*getFurthestDistance3Approx(_bidders, _items);
    maxVal = pow(maxVal, wassersteinPower);
    weightAdjConst = maxVal;
}

DebugOptimalBid AuctionOracleKDTree::getOptimalBidDebug(IdxType bidderIdx)
{
    DebugOptimalBid result;
    DiagramPoint bidder = bidders[bidderIdx];
    DnnPoint bidderDnn;
    bidderDnn[0] = bidder.getRealX();
    bidderDnn[1] = bidder.getRealY();

    std::vector<IdxValPair> candItems;
    

    if ( bidder.isDiagonal() ) {
        // 
        auto twoBestItems = kdtree->findK(bidderDnn, 2);
        candItems.push_back( std::make_pair(twoBestItems[0].p->id(), twoBestItems[0].d) );
        candItems.push_back( std::make_pair(twoBestItems[1].p->id(), twoBestItems[1].d) );
        assert(diagItemsHeap.size() > 1);
        auto topDiagIter = diagItemsHeap.ordered_begin();
        auto topDiag1 = *topDiagIter++;
        auto topDiag2 = *topDiagIter;
        candItems.push_back(topDiag1);
        candItems.push_back(topDiag2);
        assert(candItems.size() == 4);
        std::sort(candItems.begin(), candItems.end(), CompPairsBySecondStruct());
        assert(candItems[3].second >= candItems[2].second);
        assert(candItems[2].second >= candItems[1].second);
        assert(candItems[1].second >= candItems[0].second);
    } else {
        auto twoBestItems = kdtreeAll->findK(bidderDnn, 2);
        candItems.push_back( std::make_pair(twoBestItems[0].p->id(), twoBestItems[0].d) );
        candItems.push_back( std::make_pair(twoBestItems[1].p->id(), twoBestItems[1].d) );
        //size_t projItemIdx { biddersProjIndices.at(bidderIdx) };
        //assert(items[projItemIdx].projId == bidder.id);
        //double projItemValue { pow(distLp(bidder, items[projItemIdx]), wassersteinPower) + prices.at(projItemIdx) };
        //candItems.push_back( std::make_pair(projItemIdx, projItemValue) );
        assert(candItems.size() == 2);
        assert(candItems[1].second >= candItems[0].second);
    }

    result.bestItemIdx = candItems[0].first;
    result.secondBestItemIdx = candItems[1].first;
    result.bestItemValue = candItems[0].second;
    result.secondBestItemValue = candItems[1].second;
    //double bestItemsPrice = prices[bestItemIdx];
    //if (items[result.bestItemIdx].type == DiagramPoint::DIAG) {
        //double bestItemValue1 = pow( distLp(bidder, items[result.bestItemIdx]), q) + prices[result.bestItemIdx];
        //if ( fabs(result.bestItemValue - bestItemValue1) > 1e-6 ) {
            //std::cerr << "XXX: " << result.bestItemValue << " vs " << bestItemValue1 << std::endl;
            //result.bestItemValue = bestItemValue1;
        //}

    //}


    // checking code
    /*
    
    DebugOptimalBid debugMyResult(result);
    DebugOptimalBid debugNaiveResult;
    debugNaiveResult.bestItemValue = 1e20;
    debugNaiveResult.secondBestItemValue = 1e20;
    double currItemValue;
    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        //if ( bidders[bidderIdx].type == DiagramPoint::NORMAL and
                //items[itemIdx].type == DiagramPoint::DIAG and
                //bidders[bidderIdx].projId != items[itemIdx].id)
            //continue;

        currItemValue = pow(distLp(bidders[bidderIdx], items[itemIdx]), wassersteinPower) + prices[itemIdx];
        if (currItemValue < debugNaiveResult.bestItemValue) {
            debugNaiveResult.bestItemValue = currItemValue;
            debugNaiveResult.bestItemIdx  = itemIdx;
        }
    }

    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        if (itemIdx == debugNaiveResult.bestItemIdx) {
            continue;
        }
        currItemValue = pow(distLp(bidders[bidderIdx], items[itemIdx]), wassersteinPower) + prices[itemIdx];
        if (currItemValue < debugNaiveResult.secondBestItemValue) {
            debugNaiveResult.secondBestItemValue = currItemValue;
            debugNaiveResult.secondBestItemIdx = itemIdx;
        }
    }
    //std::cout << "got naive result" << std::endl;

    if ( fabs( debugMyResult.bestItemValue - debugNaiveResult.bestItemValue ) > 1e-6 or
            fabs( debugNaiveResult.secondBestItemValue - debugMyResult.secondBestItemValue) > 1e-6 ) {
        kdtreeAll->printWeights();
        std::cerr << "bidderIdx = " << bidderIdx << "; ";
        std::cerr << bidders[bidderIdx] << std::endl;
        for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
            std::cout << itemIdx << ": " << items[itemIdx] << "; price = " << prices[itemIdx] << std::endl;
        }
        std::cerr << "debugMyResult: " << debugMyResult << std::endl;
        std::cerr << "debugNaiveResult: " << debugNaiveResult << std::endl;
        //std::cerr << "twoBestItems: " << twoBestItems[0].d << " " << twoBestItems[1].d << std::endl;
        assert(false);
    }
    //std::cout << "returning" << std::endl;

    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemValue = " << bestItemValue << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << secondBestItemIdx << "; secondBestValue = " << secondBestItemValue << "; secondBestPrice = " << prices[secondBestItemIdx] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemsDist= " << (weightAdjConst - bestItemValue) << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << secondBestItemIdx << "; secondBestDist= " << (weightAdjConst - secondBestItemValue) << "; secondBestPrice = " << prices[secondBestItemIdx] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    */

    return result;
}

IdxValPair AuctionOracleKDTree::getOptimalBid(IdxType bidderIdx)
{
    IdxValPair result;
    DebugOptimalBid debugMyResult = getOptimalBidDebug(bidderIdx);
    result.first = debugMyResult.bestItemIdx;
    result.second = ( debugMyResult.secondBestItemValue - debugMyResult.bestItemValue ) + prices[debugMyResult.bestItemIdx] + epsilon;
    return result;
}
/*
a_{ij} = d_{ij} 
value_{ij} = a_{ij} + price_j
*/
void AuctionOracleKDTree::setPrice(IdxType itemIdx, double newPrice)
{
    assert(prices.size() == items.size());
    assert( 0 < diagHeapHandles.size() and diagHeapHandles.size() <= items.size());
    assert(newPrice > prices.at(itemIdx));
    prices[itemIdx] = newPrice;
    if ( items[itemIdx].isNormal() ) {
        assert(0 <= itemIdx and itemIdx < kdtreeItems.size());
        assert(0 <= kdtreeItems[itemIdx] and kdtreeItems[itemIdx] < dnnPointHandles.size());
        kdtree->increase_weight( dnnPointHandles[kdtreeItems[itemIdx]], newPrice);
        kdtreeAll->increase_weight( dnnPointHandlesAll[itemIdx], newPrice);
    } else {
        kdtreeAll->increase_weight( dnnPointHandlesAll[itemIdx], newPrice);
        assert(diagHeapHandles.size() > heapHandlesIndices.at(itemIdx));
        diagItemsHeap.decrease(diagHeapHandles[heapHandlesIndices[itemIdx]], std::make_pair(itemIdx, newPrice));
    } 
}

void AuctionOracleKDTree::adjustPrices(void)
{
}

AuctionOracleKDTree::~AuctionOracleKDTree()
{
    delete kdtree;
    delete kdtreeAll;
}

void AuctionOracleKDTree::setEpsilon(double newVal) 
{
    assert(newVal >= 0.0);
    epsilon = newVal;
}

// *****************************
// AuctionOracleRestricted
// *****************************
AuctionOracleRestricted::AuctionOracleRestricted(const std::vector<DiagramPoint>& b, 
                                         const std::vector<DiagramPoint>& g, 
                                         double _wassersteinPower,
                                         double internal_p) :
    AuctionOracleAbstract(b, g, _wassersteinPower, internal_p),
    maxVal(0.0)
{
    assert(b.size() == g.size() );
    assert(b.size() > 1);

    weightMatrix.reserve(b.size());
    for(const auto& pointA : bidders) {
        std::vector<double> weightVec;
        weightVec.clear();
        weightVec.reserve(b.size());
        for(const auto& pointB : items) {
            double val = pow(distLp(pointA, pointB, internal_p), wassersteinPower);
            if (val > maxVal) {
                maxVal = val;
            }
            weightVec.push_back( val );
        }
        weightMatrix.push_back(weightVec);
    }
}

IdxValPair AuctionOracleRestricted::getOptimalBid(const IdxType bidderIdx) 
{
    assert(bidderIdx >=0 and bidderIdx < static_cast<IdxType>(bidders.size()) );

    const auto bidder = bidders[bidderIdx];
    
    IdxType bestItemIdx { -1 };
    double bestItemValue { std::numeric_limits<double>::max() };
    //IdxType secondBestItemIdx { -1 };
    double secondBestItemValue { std::numeric_limits<double>::max() };

    // find best items idx
    for(IdxType itemIdx = 0; itemIdx < static_cast<IdxType>(items.size()); ++itemIdx) {
        // non-diagonal point should be matched either to another
        // non-diagonal point or to its own projection
        if (isRestricted and bidder.isNormal() ) {
            auto item = items[itemIdx];
            if (item.isDiagonal() and itemIdx != bidderIdx)
                continue;
        }
        auto currItemValue = weightMatrix[bidderIdx][itemIdx] + prices[itemIdx];
        if ( currItemValue < bestItemValue ) {
            bestItemValue = currItemValue;
            bestItemIdx = itemIdx;
        }
    }

    // find second best items idx and value

    for(IdxType itemIdx = 0; itemIdx < static_cast<IdxType>(items.size()); ++itemIdx) {
        // non-diagonal point should be matched either to another
        // non-diagonal point or to its own projection
        if (isRestricted and bidder.isNormal() ) {
            auto itemsItem = items[itemIdx];
            if (itemsItem.isDiagonal() and itemIdx != bidderIdx)
                continue;
        }

        if (static_cast<IdxType>(itemIdx) == bestItemIdx)
            continue;

        auto currItemValue = weightMatrix[bidderIdx][itemIdx] + prices[itemIdx];
        if ( currItemValue < secondBestItemValue ) {
            secondBestItemValue = currItemValue;
            //secondBestItemIdx = itemIdx;
        }
    }

    assert(bestItemValue <= secondBestItemValue);

    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemValue = " << bestItemValue << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestValue = " << secondBestItemValue << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemValue = " << bestItemValue << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << secondBestItemIdx << "; secondBestValue = " << secondBestItemValue << "; secondBestPrice = " << prices[secondBestItemIdx] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemsDist= " << (weightAdjConst -  bestItemValue) << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestDist= " << (weightAdjConst - secondBestItemValue) << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;

    // bid value: price + value difference + epsilon
    
    return std::make_pair(bestItemIdx, 
                          prices[bestItemIdx] + 
                          ( -bestItemValue + secondBestItemValue ) +
                          epsilon );
}

void AuctionOracleRestricted::setPrice(const IdxType itemIdx, const double newPrice)
{
    assert(prices.at(itemIdx) < newPrice );
    prices[itemIdx] = newPrice;
}

// *****************************
// AuctionOracleKDTreeRestricted
// *****************************

AuctionOracleKDTreeRestricted::AuctionOracleKDTreeRestricted(const std::vector<DiagramPoint>& _bidders, 
        const std::vector<DiagramPoint>& _items, 
        const double _wassersteinPower,
        const double internal_p) :
    AuctionOracleAbstract(_bidders, _items, _wassersteinPower, internal_p),
    heapHandlesIndices(items.size(), std::numeric_limits<size_t>::max()),
    kdtreeItems(items.size(), std::numeric_limits<size_t>::max()),
    bestDiagonalItemsComputed(false)
{
    size_t dnnItemIdx { 0 };
    size_t trueIdx { 0 };
    dnnPoints.clear();
    // store normal items in kd-tree
    for(const auto& g : items) {
        if (g.isNormal() ) {
            kdtreeItems[trueIdx] = dnnItemIdx;
            // index of items is id of dnn-point
            DnnPoint p(trueIdx);
            p[0] = g.getRealX();
            p[1] = g.getRealY();
            dnnPoints.push_back(p);
            assert(dnnItemIdx == dnnPoints.size() - 1);
            dnnItemIdx++;
        }
        trueIdx++;
    }

    assert(dnnPoints.size() < items.size() );
    for(size_t i = 0; i < dnnPoints.size(); ++i) {
        dnnPointHandles.push_back(&dnnPoints[i]);
    }
    DnnTraits traits;
    traits.internal_p = internal_p;
    kdtree = new dnn::KDTree<DnnTraits>(traits, dnnPointHandles, wassersteinPower);
    
    size_t handleIdx {0};
    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        if (items[itemIdx].isDiagonal()) {
            heapHandlesIndices[itemIdx] = handleIdx++;
            diagHeapHandles.push_back(diagItemsHeap.push(std::make_pair(itemIdx, 0)));
        }
    }
    //to-do: remove maxVal from 
    maxVal = 3*getFurthestDistance3Approx(_bidders, _items);
    maxVal = pow(maxVal, wassersteinPower);
    weightAdjConst = maxVal;
}

DebugOptimalBid AuctionOracleKDTreeRestricted::getOptimalBidDebug(IdxType bidderIdx)
{
    DebugOptimalBid result;
    DiagramPoint bidder = bidders[bidderIdx];

    // corresponding point is always considered as a candidate
    // if bidder is a diagonal point, projItem is a normal point, 
    // and vice versa.

    size_t projItemIdx = bidderIdx;
    assert( 0 <= projItemIdx and projItemIdx < items.size() );
    DiagramPoint projItem = items[projItemIdx];
    assert(projItem.type != bidder.type);
    //assert(projItem.projId == bidder.id);
    //assert(projItem.id == bidder.projId);
    // todo: store precomputed distance?
    double projItemValue = pow(distLp(bidder, projItem, internal_p), wassersteinPower) + prices[projItemIdx];
   
    if (bidder.isDiagonal()) {
        // for diagonal bidder the only normal point has already been added
        // the other 2 candidates are diagonal items only, get from the heap
        // with prices
        assert(diagItemsHeap.size() > 1);
        if (!bestDiagonalItemsComputed) {
            auto topDiagIter = diagItemsHeap.ordered_begin();
            bestDiagonalItemIdx = topDiagIter->first;
            bestDiagonalItemValue = topDiagIter->second;
            topDiagIter++;
            secondBestDiagonalItemIdx = topDiagIter->first;
            secondBestDiagonalItemValue = topDiagIter->second;
            bestDiagonalItemsComputed = true;
        }

        if ( projItemValue < bestDiagonalItemValue) {
            result.bestItemIdx = projItemIdx;
            result.bestItemValue = projItemValue;
            result.secondBestItemIdx = bestDiagonalItemIdx;
            result.secondBestItemValue = bestDiagonalItemValue;
        } else if (projItemValue < secondBestDiagonalItemValue) {
            result.bestItemIdx = bestDiagonalItemIdx;
            result.bestItemValue = bestDiagonalItemValue;
            result.secondBestItemIdx = projItemIdx;
            result.secondBestItemValue = projItemValue;
        } else {
            result.bestItemIdx = bestDiagonalItemIdx;
            result.bestItemValue = bestDiagonalItemValue;
            result.secondBestItemIdx = secondBestDiagonalItemIdx;
            result.secondBestItemValue = secondBestDiagonalItemValue;
        }
    } else {
        // for normal bidder get 2 best items among non-diagonal points from
        // kdtree
        DnnPoint bidderDnn;
        bidderDnn[0] = bidder.getRealX();
        bidderDnn[1] = bidder.getRealY();
        auto twoBestItems = kdtree->findK(bidderDnn, 2);
        size_t bestNormalItemIdx { twoBestItems[0].p->id() };
        double bestNormalItemValue { twoBestItems[0].d };
        size_t secondBestNormalItemIdx { twoBestItems[1].p->id() };
        double secondBestNormalItemValue { twoBestItems[1].d };

        if ( projItemValue < bestNormalItemValue) {
            result.bestItemIdx = projItemIdx;
            result.bestItemValue = projItemValue;
            result.secondBestItemIdx = bestNormalItemIdx;
            result.secondBestItemValue = bestNormalItemValue;
        } else if (projItemValue < secondBestNormalItemValue) {
            result.bestItemIdx = bestNormalItemIdx;
            result.bestItemValue = bestNormalItemValue;
            result.secondBestItemIdx = projItemIdx;
            result.secondBestItemValue = projItemValue;
        } else {
            result.bestItemIdx = bestNormalItemIdx;
            result.bestItemValue = bestNormalItemValue;
            result.secondBestItemIdx = secondBestNormalItemIdx;
            result.secondBestItemValue = secondBestNormalItemValue;
        }
    }

    return result;

    //std::cout << "got result: " << result << std::endl;
    //double bestItemsPrice = prices[bestItemIdx];
    //if (items[result.bestItemIdx].type == DiagramPoint::DIAG) {
        //double bestItemValue1 = pow( distLp(bidder, items[result.bestItemIdx]), wassersteinPower) + prices[result.bestItemIdx];
        //if ( fabs(result.bestItemValue - bestItemValue1) > 1e-6 ) {
            //std::cerr << "XXX: " << result.bestItemValue << " vs " << bestItemValue1 << std::endl;
            //result.bestItemValue = bestItemValue1;
        //}

    //}


    // checking code
    
    /*
    DebugOptimalBid debugMyResult(result);
    DebugOptimalBid debugNaiveResult;
    debugNaiveResult.bestItemValue = 1e20;
    debugNaiveResult.secondBestItemValue = 1e20;
    double currItemValue;
    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        //if ( bidders[bidderIdx].type == DiagramPoint::NORMAL and
                //items[itemIdx].type == DiagramPoint::DIAG and
                //bidders[bidderIdx].projId != items[itemIdx].id)
            //continue;

        currItemValue = pow(distLp(bidders[bidderIdx], items[itemIdx]), wassersteinPower) + prices[itemIdx];
        if (currItemValue < debugNaiveResult.bestItemValue) {
            debugNaiveResult.bestItemValue = currItemValue;
            debugNaiveResult.bestItemIdx  = itemIdx;
        }
    }

    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        if (itemIdx == debugNaiveResult.bestItemIdx) {
            continue;
        }
        currItemValue = pow(distLp(bidders[bidderIdx], items[itemIdx]), wassersteinPower) + prices[itemIdx];
        if (currItemValue < debugNaiveResult.secondBestItemValue) {
            debugNaiveResult.secondBestItemValue = currItemValue;
            debugNaiveResult.secondBestItemIdx = itemIdx;
        }
    }
    //std::cout << "got naive result" << std::endl;

    if ( fabs( debugMyResult.bestItemValue - debugNaiveResult.bestItemValue ) > 1e-6 or
            fabs( debugNaiveResult.secondBestItemValue - debugMyResult.secondBestItemValue) > 1e-6 ) {
        kdtreeAll->printWeights();
        std::cerr << "bidderIdx = " << bidderIdx << "; ";
        std::cerr << bidders[bidderIdx] << std::endl;
        for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
            std::cout << itemIdx << ": " << items[itemIdx] << "; price = " << prices[itemIdx] << std::endl;
        }
        std::cerr << "debugMyResult: " << debugMyResult << std::endl;
        std::cerr << "debugNaiveResult: " << debugNaiveResult << std::endl;
        //std::cerr << "twoBestItems: " << twoBestItems[0].d << " " << twoBestItems[1].d << std::endl;
        assert(false);
    }
    //std::cout << "returning" << std::endl;

    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemValue = " << bestItemValue << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << secondBestItemIdx << "; secondBestValue = " << secondBestItemValue << "; secondBestPrice = " << prices[secondBestItemIdx] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemsDist= " << (weightAdjConst - bestItemValue) << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << secondBestItemIdx << "; secondBestDist= " << (weightAdjConst - secondBestItemValue) << "; secondBestPrice = " << prices[secondBestItemIdx] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    */
    return result;
}

IdxValPair AuctionOracleKDTreeRestricted::getOptimalBid(IdxType bidderIdx)
{

    
    DiagramPoint bidder = bidders[bidderIdx];

    // corresponding point is always considered as a candidate
    // if bidder is a diagonal point, projItem is a normal point, 
    // and vice versa.
    
    size_t bestItemIdx;
    double bestItemValue;
    double secondBestItemValue;


    size_t projItemIdx = bidderIdx;
    assert( 0 <= projItemIdx and projItemIdx < items.size() );
    DiagramPoint projItem = items[projItemIdx];
    assert(projItem.type != bidder.type);
    //assert(projItem.projId == bidder.id);
    //assert(projItem.id == bidder.projId);
    // todo: store precomputed distance?
    double projItemValue = pow(distLp(bidder, projItem, internal_p), wassersteinPower) + prices[projItemIdx];
   
    if (bidder.isDiagonal()) {
        // for diagonal bidder the only normal point has already been added
        // the other 2 candidates are diagonal items only, get from the heap
        // with prices

        if (not bestDiagonalItemsComputed) {
            auto topDiagIter = diagItemsHeap.ordered_begin();
            bestDiagonalItemIdx = topDiagIter->first;
            bestDiagonalItemValue = topDiagIter->second;
            if (diagItemsHeap.size() > 1) {
                topDiagIter++;
                secondBestDiagonalItemIdx = topDiagIter->first;
                secondBestDiagonalItemValue = topDiagIter->second;
            } else {
                // there is only one diagonal point at all,
                // ensure that second best diagonal value
                // will lose to projection item
                secondBestDiagonalItemValue = std::numeric_limits<double>::max();
                secondBestDiagonalItemIdx = -1;
            }
            bestDiagonalItemsComputed = true;
        }

        if ( projItemValue < bestDiagonalItemValue) {
            bestItemIdx = projItemIdx;
            bestItemValue = projItemValue;
            secondBestItemValue = bestDiagonalItemValue;
        } else if (projItemValue < secondBestDiagonalItemValue) {
            bestItemIdx = bestDiagonalItemIdx;
            bestItemValue = bestDiagonalItemValue;
            secondBestItemValue = projItemValue;
        } else {
            bestItemIdx = bestDiagonalItemIdx;
            bestItemValue = bestDiagonalItemValue;
            secondBestItemValue = secondBestDiagonalItemValue;
        }
    } else {
        // for normal bidder get 2 best items among non-diagonal points from
        // kdtree
        DnnPoint bidderDnn;
        bidderDnn[0] = bidder.getRealX();
        bidderDnn[1] = bidder.getRealY();
        auto twoBestItems = kdtree->findK(bidderDnn, 2);
        size_t bestNormalItemIdx { twoBestItems[0].p->id() };
        double bestNormalItemValue { twoBestItems[0].d };
        // if there is only one off-diagonal point in the second diagram,
        // kd-tree will not return the second candidate. 
        // Set its value to inf, so it will always lose to the value of the projection
        double secondBestNormalItemValue { twoBestItems.size() == 1 ? std::numeric_limits<double>::max() : twoBestItems[1].d };

        if ( projItemValue < bestNormalItemValue) {
            bestItemIdx = projItemIdx;
            bestItemValue = projItemValue;
            secondBestItemValue = bestNormalItemValue;
        } else if (projItemValue < secondBestNormalItemValue) {
            bestItemIdx = bestNormalItemIdx;
            bestItemValue = bestNormalItemValue;
            secondBestItemValue = projItemValue;
        } else {
            bestItemIdx = bestNormalItemIdx;
            bestItemValue = bestNormalItemValue;
            secondBestItemValue = secondBestNormalItemValue;
        }
    }

    IdxValPair result;

    assert( secondBestItemValue >= bestItemValue );

    result.first = bestItemIdx;
    result.second = ( secondBestItemValue - bestItemValue ) + prices[bestItemIdx] + epsilon;
    return result;
}
/*
a_{ij} = d_{ij} 
value_{ij} = a_{ij} + price_j
*/
void AuctionOracleKDTreeRestricted::setPrice(IdxType itemIdx, double newPrice)
{
    assert(prices.size() == items.size());
    assert( 0 < diagHeapHandles.size() and diagHeapHandles.size() <= items.size());
    assert(newPrice > prices.at(itemIdx));
    prices[itemIdx] = newPrice;
    if ( items[itemIdx].isNormal() ) {
        assert(0 <= itemIdx and itemIdx < kdtreeItems.size());
        assert(0 <= kdtreeItems[itemIdx] and kdtreeItems[itemIdx] < dnnPointHandles.size());
        kdtree->increase_weight( dnnPointHandles[kdtreeItems[itemIdx]], newPrice);
    } else {
        assert(diagHeapHandles.size() > heapHandlesIndices.at(itemIdx));
        diagItemsHeap.decrease(diagHeapHandles[heapHandlesIndices[itemIdx]], std::make_pair(itemIdx, newPrice));
        bestDiagonalItemsComputed = false;
    }
}

void AuctionOracleKDTreeRestricted::adjustPrices(void)
{
}

AuctionOracleKDTreeRestricted::~AuctionOracleKDTreeRestricted()
{
    delete kdtree;
}

void AuctionOracleKDTreeRestricted::setEpsilon(double newVal) 
{
    assert(newVal >= 0.0);
    epsilon = newVal;
}

std::ostream& operator<< (std::ostream& output, const DebugOptimalBid& db)
{
    output << "bestItemValue = " << db.bestItemValue;
    output << "; bestItemIdx = " << db.bestItemIdx;
    output << "; secondBestItemValue = " << db.secondBestItemValue;
    output << "; secondBestItemIdx = " << db.secondBestItemIdx;
    return output;
}

} // end of namespace geom_ws
