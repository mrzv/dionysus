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

#ifndef HERA_BASIC_DEFS_BT_H
#define HERA_BASIC_DEFS_BT_H

#ifdef _WIN32
#include <ciso646>
#endif

#include <utility>
#include <vector>
#include <stdexcept>
#include <math.h>
#include <cstddef>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <assert.h>

#include "def_debug_bt.h"
#include <hera/common.h>
#include <iostream>

namespace hera {

    namespace bt {

        typedef int IdType;

//        template<class Real = double>
//        struct DiagramPoint
//        {
//            // Points above the diagonal have type NORMAL
//            // Projections onto the diagonal have type DIAG
//            // for DIAG points only x-coordinate is relevant
//            // to-do: add getters/setters, checks in constructors, etc
//            enum Type
//            {
//                NORMAL, DIAG
//            };
//            // data members
//        private:
//            Real x, y;
//        public:
//            Type type;
//            IdType id;
//            IdType user_id;

//            // operators, constructors
//            bool operator==(const DiagramPoint<Real>& other) const
//            {
//                // compare by id only
//                assert(this->id >= MinValidId);
//                assert(other.id >= MinValidId);
//                bool areEqual { this->id == other.id };
//                assert(!areEqual or ((this->x == other.x) and (this->y == other.y) and (this->type == other.type)));
//                return areEqual;
//            }

//            bool operator!=(const DiagramPoint& other) const
//            {
//                return !(*this == other);
//            }

//            DiagramPoint() :
//                    x(0.0),
//                    y(0.0),
//                    type(DiagramPoint::DIAG),
//                    id(MinValidId - 1),
//                    user_id(-1)
//            {
//            }

//            DiagramPoint(Real _x, Real _y, Type _type, IdType _id, IdType _user_id) :
//                    x(_x),
//                    y(_y),
//                    type(_type),
//                    id(_id),
//                    user_id(_user_id)
//            {
//                if (_y == _x and _type != DIAG) {
//                    throw std::runtime_error("Point on the main diagonal must have DIAG type");
//                }

//            }


//            bool is_diagonal() const
//            { return type == DIAG; }

//            bool is_normal() const
//            { return type == NORMAL; }

//            bool is_infinity() const
//            {
//                return x == std::numeric_limits<Real>::infinity() or
//                       x == -std::numeric_limits<Real>::infinity() or
//                       y == std::numeric_limits<Real>::infinity() or
//                       y == -std::numeric_limits<Real>::infinity();
//            }

//            Real inline getRealX() const // return the x-coord
//            {
//                return x;
//            }

//            Real inline getRealY() const // return the y-coord
//            {
//                return y;
//            }

//            IdType inline get_user_id() const
//            {
//                if (is_normal())
//                    return user_id;
//                else
//                    return -1;
//            }

//            Real inline get_persistence(const Real internal_p = get_infinity()) const
//            {
//                if (is_diagonal())
//                    return 0.0;
//                Real pers = fabs(y - x) / 2;
//                if (internal_p == get_infinity()) {
//                    return pers;
//                } else if (internal_p == 1.0) {
//                    return 2 * pers;
//                } else {
//                    return std::pow(static_cast<Real>(2), static_cast<Real>(1) / internal_p);
//                }
//            }

//#ifndef FOR_R_TDA

//            friend std::ostream& operator<<(std::ostream& output, const DiagramPoint& p)
//            {
//                if (p.is_diagonal()) {
//                    output << "(" << p.x << ", " << p.y << ", " << 0.5 * (p.x + p.y) << ", " << p.id << " DIAG )";
//                } else {
//                    output << "(" << p.x << ", " << p.y << ", " << p.id << " NORMAL)";
//                }
//                return output;
//            }
//#endif
//        };

        template<class Real>
        using MatchingEdge = std::pair<DiagramPoint<Real>, DiagramPoint<Real>>;

        // this function works with points at infinity as well
        // not needed in actual computation, since these points are processed
        // separately, but is useful in tests
        template<class Real>
        inline Real dist_l_inf_slow(const DiagramPoint<Real>& a, const DiagramPoint<Real>& b)
        {
            if (a.is_diagonal() and b.is_diagonal()) {
                // distance between points on the diagonal is 0
                return 0.0;
            }
            // otherwise distance is a usual l-inf distance
            //
            Real dx = (a.getRealX() == b.getRealX()) ? 0.0 : fabs(a.getRealX() - b.getRealX());
            Real dy = (a.getRealY() == b.getRealY()) ? 0.0 : fabs(a.getRealY() - b.getRealY());
            Real result = std::max(dx, dy);
            if (std::isnan(result))
                result = std::numeric_limits<Real>::infinity();
            return result;
        }


        //template<class Real = double>
        //typedef std::unordered_set<Point, PointHash> PointSet;
        template<class Real_ = double>
        class DiagramPointSet;

        template<class Real>
        void addProjections(DiagramPointSet<Real>& A, DiagramPointSet<Real>& B);

        template<class Real_>
        class DiagramPointSet
        {
        public:

            using Real = Real_;
            using DgmPoint = DiagramPoint<Real>;
            using DgmPointHash = DiagramPointHash<Real>;
            using const_iterator = typename std::unordered_set<DgmPoint, DgmPointHash>::const_iterator;
            using iterator = typename std::unordered_set<DgmPoint, DgmPointHash>::iterator;

        private:

            bool isLinked { false };
            IdType maxId { 1 };
            std::unordered_set<DgmPoint, DgmPointHash> points;

        public:

            void insert(const DgmPoint& p)
            {
                points.insert(p);
                if (p.id > maxId) {
                    maxId = p.id + 1;
                }
            }

            void erase(const DgmPoint& p, bool doCheck = true)
            {
                // if doCheck, erasing non-existing elements causes assert
                auto it = points.find(p);
                if (it != points.end()) {
                    points.erase(it);
                } else {
                    assert(!doCheck);
                }
                (void)doCheck; // to suppress warning, doCheck is unused in Release mode
            }


            void erase(const const_iterator it)
            {
                points.erase(it);
            }

            void removeDiagonalPoints()
            {
                if (isLinked) {
                    auto ptIter = points.begin();
                    while (ptIter != points.end()) {
                        if (ptIter->is_diagonal()) {
                            ptIter = points.erase(ptIter);
                        } else {
                            ptIter++;
                        }
                    }
                    isLinked = false;
                }
            }

            size_t size() const
            {
                return points.size();
            }

            void reserve(const size_t newSize)
            {
                points.reserve(newSize);
            }

            void clear()
            {
                points.clear();
            }

            bool empty() const
            {
                return points.empty();
            }

            bool hasElement(const DgmPoint& p) const
            {
                return points.find(p) != points.end();
            }

            iterator find(const DgmPoint& p)
            {
                return points.find(p);
            }

            iterator begin()
            {
                return points.begin();
            }

            iterator end()
            {
                return points.end();
            }

            const_iterator cbegin() const
            {
                return points.cbegin();
            }

            const_iterator cend() const
            {
                return points.cend();
            }


            const_iterator find(const DgmPoint& p) const
            {
                return points.find(p);
            }


            friend std::ostream& operator<<(std::ostream& output, const DiagramPointSet& ps)
            {
                output << "{ ";
                for (auto pit = ps.cbegin(); pit != ps.cend(); ++pit) {
                    output << *pit << ", ";
                }
                output << "\b\b }";
                return output;
            }


            friend void addProjections<Real_>(DiagramPointSet<Real_>& A, DiagramPointSet<Real_>& B);

            template<class PairIterator>
            void fillIn(PairIterator begin_iter, PairIterator end_iter)
            {
                isLinked = false;
                clear();
                IdType uniqueId = 0;
                IdType user_id = 0;
                for (auto iter = begin_iter; iter != end_iter; ++iter) {
                    insert(DgmPoint(iter->first, iter->second, DgmPoint::NORMAL, uniqueId++, user_id++));
                }
            }

            template<class PointContainer>
            void fillIn(const PointContainer& dgm_cont)
            {
                using Traits = DiagramTraits<PointContainer>;
                isLinked = false;
                clear();
                IdType uniqueId = 0;
                IdType user_id = 0;
                for (const auto& pt : dgm_cont) {
                    Real x = Traits::get_x(pt);
                    Real y = Traits::get_y(pt);
                    insert(DgmPoint(x, y, DgmPoint::NORMAL, uniqueId++, user_id++));
                }
            }


            // ctor from range
            template<class PairIterator>
            DiagramPointSet(PairIterator begin_iter, PairIterator end_iter)
            {
                fillIn(begin_iter, end_iter);
            }

            // ctor from container, uses DiagramTraits
            template<class PointContainer>
            DiagramPointSet(const PointContainer& dgm)
            {
                fillIn(dgm);
            }


            // default ctor, empty diagram
            DiagramPointSet(IdType minId = 1) :
                    maxId(minId + 1)
            {};

            IdType nextId()
            { return maxId + 1; }

        }; // DiagramPointSet


        template<class Real, class DiagPointContainer>
        Real getFurthestDistance3Approx(DiagPointContainer& A, DiagPointContainer& B)
        {
            Real result { 0.0 };
            DiagramPoint<Real> begA = *(A.begin());
            DiagramPoint<Real> optB = *(B.begin());
            for (const auto& pointB : B) {
                if (dist_l_inf(begA, pointB) > result) {
                    result = dist_l_inf(begA, pointB);
                    optB = pointB;
                }
            }
            for (const auto& pointA : A) {
                if (dist_l_inf(pointA, optB) > result) {
                    result = dist_l_inf(pointA, optB);
                }
            }
            return result;
        }

        // preprocess diagrams A and B by adding projections onto diagonal of points of
        // A to B and vice versa. Also removes points at infinity!
        // NB: ids of points will be changed!
        template<class Real_>
        void addProjections(DiagramPointSet<Real_>& A, DiagramPointSet<Real_>& B)
        {

            using Real = Real_;
            using DgmPoint = DiagramPoint<Real>;
            using DgmPointSet = DiagramPointSet<Real>;

            IdType uniqueId { 0 };
            DgmPointSet newA, newB;

            // copy normal points from A to newA
            // add projections to newB
            for (auto& pA : A) {
                if (pA.is_normal() and not pA.is_infinity()) {
                    // add pA's projection to B
                    DgmPoint dpA { pA.getRealX(), pA.getRealY(), DgmPoint::NORMAL, uniqueId++, pA.user_tag };
                    DgmPoint dpB { (pA.getRealX() + pA.getRealY()) / 2, (pA.getRealX() + pA.getRealY()) / 2,
                                   DgmPoint::DIAG, uniqueId++, -1 - pA.user_tag };
                    newA.insert(dpA);
                    newB.insert(dpB);
                }
            }

            for (auto& pB : B) {
                if (pB.is_normal() and not pB.is_infinity()) {
                    // add pB's projection to A
                    DgmPoint dpB { pB.getRealX(), pB.getRealY(), DgmPoint::NORMAL, uniqueId++, pB.user_tag };
                    DgmPoint dpA { (pB.getRealX() + pB.getRealY()) / 2, (pB.getRealX() + pB.getRealY()) / 2,
                                   DgmPoint::DIAG, uniqueId++, -1 - pB.user_tag };
                    newB.insert(dpB);
                    newA.insert(dpA);
                }
            }

            A = newA;
            B = newB;
            A.isLinked = true;
            B.isLinked = true;
        }

        //#ifndef FOR_R_TDA

        //template<class Real>
        //std::ostream& operator<<(std::ostream& output, const DiagramPoint<Real>& p)
        //{
        //    if ( p.is_diagonal() ) {
        //        output << "(" << p.x << ", " << p.y << ", " <<  0.5 * (p.x + p.y) << ", "  << p.id << " DIAG )";
        //    } else {
        //        output << "(" << p.x << ", " << p.y << ", " << p.id << " NORMAL)";
        //    }
        //    return output;
        //}

        //template<class Real>
        //std::ostream& operator<<(std::ostream& output, const DiagramPointSet<Real>& ps)
        //{
        //    output << "{ ";
        //    for(auto pit = ps.cbegin(); pit != ps.cend(); ++pit) {
        //        output << *pit << ", ";
        //    }
        //    output << "\b\b }";
        //    return output;
        //}
        //#endif // FOR_R_TDA


    } // end namespace bt
} // end namespace hera
#endif  // HERA_BASIC_DEFS_BT_H
