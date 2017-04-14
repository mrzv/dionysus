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

#ifndef AUCTION_ORACLE_H
#define AUCTION_ORACLE_H


#define USE_BOOST_HEAP

#include <map>
#include <memory>
#include <set>
#include <list>

#ifdef USE_BOOST_HEAP
#include <boost/heap/d_ary_heap.hpp>
#endif

#include "basic_defs_ws.h"
#include "dnn/geometry/euclidean-fixed.h"
#include "dnn/local/kd-tree.h"

namespace geom_ws {

struct CompPairsBySecondStruct {
    bool operator()(const IdxValPair& a, const IdxValPair& b) const
    {
        return a.second < b.second;
    }
};

// 
struct CompPairsBySecondGreaterStruct {
    bool operator()(const IdxValPair& a, const IdxValPair& b) const
    {
        return a.second > b.second;
    }
};

struct CompPairsBySecondLexStruct {
    bool operator()(const IdxValPair& a, const IdxValPair& b) const
    {
        return a.second < b.second or (a.second == b.second and a.first > b.first);
    }
};

struct CompPairsBySecondLexGreaterStruct {
    bool operator()(const IdxValPair& a, const IdxValPair& b) const
    {
        return a.second > b.second or (a.second == b.second and a.first > b.first);
    }
};

using ItemsTimePair = std::pair<IdxType, int>;

using UpdateList = std::list<ItemsTimePair>;
using UpdateListIter = UpdateList::iterator;


#ifdef USE_BOOST_HEAP
using LossesHeap = boost::heap::d_ary_heap<IdxValPair, boost::heap::arity<2>, boost::heap::mutable_<true>, boost::heap::compare<CompPairsBySecondGreaterStruct>>;
#else
template<class ComparisonStruct>
class IdxValHeap {
public:
    using InternalKeeper = std::set<IdxValPair, ComparisonStruct>;
    using handle_type = typename InternalKeeper::iterator;
    // methods
    handle_type push(const IdxValPair& val)
    {
        auto resPair = _heap.insert(val);
        assert(resPair.second);
        assert(resPair.first != _heap.end());
        return resPair.first;
    }

    void decrease(handle_type& handle, const IdxValPair& newVal)
    {
        _heap.erase(handle);
        handle = push(newVal);
    }

    size_t size() const 
    { 
        return _heap.size();
    }

    handle_type ordered_begin() 
    {
        return _heap.begin();
    }

    handle_type ordered_end() 
    {
        return _heap.end();
    }


private:
    std::set<IdxValPair, ComparisonStruct> _heap;
};

// if we store losses, the minimal value should come first
using LossesHeap = IdxValHeap<CompPairsBySecondLexStruct>;
#endif

struct DebugOptimalBid {
    DebugOptimalBid() : bestItemIdx(-1), bestItemValue(-666.666), secondBestItemIdx(-1), secondBestItemValue(-666.666) {};
    IdxType bestItemIdx;
    double bestItemValue;
    IdxType secondBestItemIdx;
    double secondBestItemValue;
};

struct AuctionOracleAbstract {
    AuctionOracleAbstract(const std::vector<DiagramPoint>& _bidders, const std::vector<DiagramPoint>& _items, const double _wassersteinPower, const double _internal_p = std::numeric_limits<double>::infinity());
    ~AuctionOracleAbstract() {}
    virtual IdxValPair getOptimalBid(const IdxType bidderIdx) = 0;
    virtual void setPrice(const IdxType itemsIdx, const double newPrice) = 0;
    virtual void adjustPrices(void) = 0;
    double getEpsilon() { return epsilon; };
    virtual void setEpsilon(double newEpsilon) { assert(newEpsilon >= 0.0); epsilon = newEpsilon; };
    std::vector<double> getPrices() { return prices; }
protected:
    const std::vector<DiagramPoint>& bidders;
    const std::vector<DiagramPoint>& items;
    std::vector<double> prices;
    double wassersteinPower;
    double epsilon;
    double internal_p;
    double getValueForBidder(size_t bidderIdx, size_t itemsIdx);
};

struct AuctionOracleLazyHeap final : AuctionOracleAbstract {
    AuctionOracleLazyHeap(const std::vector<DiagramPoint>& bidders, const std::vector<DiagramPoint>& items, const double wassersteinPower, const double _internal_p = std::numeric_limits<double>::infinity());
    ~AuctionOracleLazyHeap();
    // data members
    // temporarily make everything public
    std::vector<std::vector<double>> weightMatrix;
    //double weightAdjConst;
    double maxVal;
    // vector of heaps to find the best items
    std::vector<LossesHeap*> lossesHeap;
    std::vector<std::vector<LossesHeap::handle_type>> lossesHeapHandles;
    // methods
    void fillInLossesHeap(void);
    void setPrice(const IdxType itemsIdx, const double newPrice) override final;
    IdxValPair getOptimalBid(const IdxType bidderIdx) override final;
    double getMatchingWeight(const std::vector<IdxType>& biddersToItems) const;
    void adjustPrices(void) override final;
    // to update the queue in lazy fashion
    std::vector<UpdateListIter> itemsIterators;
    UpdateList updateList;
    std::vector<int> biddersUpdateMoments;
    int updateCounter;
    void updateQueueForBidder(const IdxType bidderIdx);
    // debug
    DebugOptimalBid getOptimalBidDebug(const IdxType bidderIdx);
};

struct AuctionOracleLazyHeapRestricted final : AuctionOracleAbstract {
     AuctionOracleLazyHeapRestricted(const std::vector<DiagramPoint>& bidders, const std::vector<DiagramPoint>& items, const double wassersteinPower, const double _internal_p = std::numeric_limits<double>::infinity());
    ~AuctionOracleLazyHeapRestricted();
    // data members
    // temporarily make everything public
    std::vector<std::vector<double>> weightMatrix;
    //double weightAdjConst;
    double maxVal;
    // vector of heaps to find the best items
    std::vector<LossesHeap*> lossesHeap;
    std::vector<std::vector<size_t>> itemsIndicesForHeapHandles;
    std::vector<std::vector<LossesHeap::handle_type>> lossesHeapHandles;
    // methods
    void fillInLossesHeap(void);
    void setPrice(const IdxType itemsIdx, const double newPrice) override final;
    IdxValPair getOptimalBid(const IdxType bidderIdx) override final;
    double getMatchingWeight(const std::vector<IdxType>& biddersToItems) const;
    void adjustPrices(void) override final;
    // to update the queue in lazy fashion
    std::vector<UpdateListIter> itemsIterators;
    UpdateList updateList;
    std::vector<int> biddersUpdateMoments;
    int updateCounter;
    void updateQueueForBidder(const IdxType bidderIdx);
    LossesHeap diagItemsHeap;
    std::vector<LossesHeap::handle_type> diagHeapHandles;
    std::vector<size_t> heapHandlesIndices;
    // debug
   
    DebugOptimalBid getOptimalBidDebug(const IdxType bidderIdx);
    
    // for diagonal points
    bool bestDiagonalItemsComputed;
    size_t bestDiagonalItemIdx;
    double bestDiagonalItemValue;
    size_t secondBestDiagonalItemIdx;
    double secondBestDiagonalItemValue;
};

struct AuctionOracleKDTree final : AuctionOracleAbstract {
    typedef dnn::Point<2, double> DnnPoint;
    typedef dnn::PointTraits<DnnPoint> DnnTraits;

    AuctionOracleKDTree(const std::vector<DiagramPoint>& bidders, const std::vector<DiagramPoint>& items, const double wassersteinPower, const double _internal_p = std::numeric_limits<double>::infinity());
    ~AuctionOracleKDTree();
    // data members
    // temporarily make everything public
    double maxVal;
    double weightAdjConst;
    dnn::KDTree<DnnTraits>* kdtree;
    std::vector<DnnPoint> dnnPoints;
    std::vector<DnnPoint*> dnnPointHandles;
    dnn::KDTree<DnnTraits>* kdtreeAll;
    std::vector<DnnPoint> dnnPointsAll;
    std::vector<DnnPoint*> dnnPointHandlesAll;
    LossesHeap diagItemsHeap;
    std::vector<LossesHeap::handle_type> diagHeapHandles;
    std::vector<size_t> heapHandlesIndices;
    std::vector<size_t> kdtreeItems;
    // vector of heaps to find the best items
    void setPrice(const IdxType itemsIdx, const double newPrice) override final;
    IdxValPair getOptimalBid(const IdxType bidderIdx) override final;
    void adjustPrices(void) override final;
    // debug routines
    DebugOptimalBid getOptimalBidDebug(IdxType bidderIdx);
    void setEpsilon(double newVal) override final;
};

struct AuctionOracleKDTreeRestricted final : AuctionOracleAbstract {
    typedef dnn::Point<2, double> DnnPoint;
    typedef dnn::PointTraits<DnnPoint> DnnTraits;

    AuctionOracleKDTreeRestricted(const std::vector<DiagramPoint>& bidders, const std::vector<DiagramPoint>& items, const double wassersteinPower, const double _internal_p = std::numeric_limits<double>::infinity());
    ~AuctionOracleKDTreeRestricted();
    // data members
    // temporarily make everything public
    double maxVal;
    double weightAdjConst;
    dnn::KDTree<DnnTraits>* kdtree;
    std::vector<DnnPoint> dnnPoints;
    std::vector<DnnPoint*> dnnPointHandles;
    std::vector<DnnPoint> dnnPointsAll;
    std::vector<DnnPoint*> dnnPointHandlesAll;
    LossesHeap diagItemsHeap;
    std::vector<LossesHeap::handle_type> diagHeapHandles;
    std::vector<size_t> heapHandlesIndices;
    std::vector<size_t> kdtreeItems;
    // vector of heaps to find the best items
    void setPrice(const IdxType itemsIdx, const double newPrice) override final;
    IdxValPair getOptimalBid(const IdxType bidderIdx) override final;
    void adjustPrices(void) override final;
    // debug routines
    DebugOptimalBid getOptimalBidDebug(IdxType bidderIdx);
    void setEpsilon(double newVal) override final;


    bool bestDiagonalItemsComputed;
    size_t bestDiagonalItemIdx;
    double bestDiagonalItemValue;
    size_t secondBestDiagonalItemIdx;
    double secondBestDiagonalItemValue;
};

struct AuctionOracleRestricted final : AuctionOracleAbstract {
    AuctionOracleRestricted(const std::vector<DiagramPoint>& bidders, const std::vector<DiagramPoint>& items, const double wassersteinPower, const double _internal_p = std::numeric_limits<double>::infinity());
    IdxValPair getOptimalBid(const IdxType bidderIdx) override;
    void setPrice(const IdxType itemsIdx, const double newPrice) override;
    void adjustPrices(void) override {};
    void setEpsilon(double newEpsilon) override { assert(newEpsilon >= 0.0); epsilon = newEpsilon; };
    // data 
    std::vector<std::vector<double>> weightMatrix;
    double maxVal;
    constexpr static bool isRestricted = true;
};

std::ostream& operator<< (std::ostream& output, const DebugOptimalBid& db);

} // end of namespace geom_ws

#endif
