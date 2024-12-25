/*

Copyright (c) 2022, M. Kerber, D. Morozov, A. Nigmetov
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

#ifndef HERA_DIAGRAM_POINT_H
#define HERA_DIAGRAM_POINT_H

#include <ostream>
#include <cassert>
#include <limits>
#include <functional>
#include <tuple>

#include "infinity.h"


namespace hera {

//    enum class OwnerType { k_none, k_normal, k_diagonal };
//
//    inline std::ostream& operator<<(std::ostream& s, const OwnerType t)
//    {
//        switch(t)
//        {
//        case OwnerType::k_none : s << "NONE"; break;
//        case OwnerType::k_normal: s << "NORMAL"; break;
//        case OwnerType::k_diagonal: s << "DIAGONAL"; break;
//        }
//        return s;
//    }

    template<class Real_ = double>
    struct DiagramPoint {
        // types
        using Real = Real_;

        // Points above the diagonal have type NORMAL
        // Projections onto the diagonal have type DIAG
        // for DIAG points only x-coordinate is relevant
        enum Type { NORMAL, DIAG };

        // data members

        Real x, y;
        Type type;
        int id;
        int user_tag;

        // methods
        DiagramPoint(Real x, Real y, Type type, int id, int user_tag = 0)
                :
                x(x), y(y), type(type), id(id), user_tag(user_tag)
        {
        };

        DiagramPoint(Real x, Real y, int id, int user_tag = 0)
                :
                x(x), y(y), type( (x == y) ? DIAG : NORMAL), id(id), user_tag(user_tag)
        {
        };

        DiagramPoint()
                :x(0), y(0), type(DIAG), id(0), user_tag(0) { };

        bool is_diagonal() const { return type == DIAG; }
        bool is_normal() const { return type == NORMAL; }

        int get_id() const { return id; }

        // TODO: here we assume that Real has true infinity
        bool is_infinity() const
        {
            return x == std::numeric_limits<Real>::infinity() or
                    x == -std::numeric_limits<Real>::infinity() or
                    y == std::numeric_limits<Real>::infinity() or
                    y == -std::numeric_limits<Real>::infinity();
        }

        Real getRealX() const // return the x-coord
        {
            if (is_normal())
                return x;
            else
                return (x + y) / 2;
        }

        Real getRealY() const // return the y-coord
        {
            if (is_normal())
                return y;
            else
                return (x + y) / 2;
        }

        Real persistence_lp(const Real p) const
        {
            if (is_diagonal())
                return 0.0;
            else {
                Real u {(getRealY() + getRealX()) / 2};
                int dim = 2;
                DiagramPoint<Real> a_proj(u, u, DiagramPoint<Real>::DIAG, /*dummy id*/ 0);
                return dist_lp(*this, a_proj, p, dim);
            }
        }

        const Real& operator[](const int idx) const
        {
            switch(idx) {
            case 0 : return x;
                break;
            case 1 : return y;
                break;
            default: throw std::out_of_range("DiagramPoint has dimension 2");
            }
        }

        Real& operator[](const int idx)
        {
            switch(idx) {
            case 0 : return x;
                break;
            case 1 : return y;
                break;
            default: throw std::out_of_range("DiagramPoint has dimension 2");
            }
        }

        // operators, constructors
        bool operator==(const DiagramPoint<Real>& other) const
        {
            return x == other.x and y == other.y and type == other.type and id == other.id and user_tag == other.user_tag;
        }

        bool operator!=(const DiagramPoint& other) const
        {
            return !(*this == other);
        }

        bool operator<(const DiagramPoint<Real>& other) const
        {
            return std::tie(type, x, y, id, user_tag) < std::tie(other.type, other.x, other.y, other.id, other.user_tag);
        }

        bool operator<=(const DiagramPoint<Real>& other) const
        {
            return std::tie(type, x, y, id, user_tag) <= std::tie(other.type, other.x, other.y, other.id, other.user_tag);
        }

        bool operator>(const DiagramPoint<Real>& other) const
        {
            return !(*this <= other);
        }

        bool operator>=(const DiagramPoint<Real>& other) const
        {
            return !(*this < other);
        }
    };

    template<class Real>
    inline std::ostream& operator<<(std::ostream& output, const DiagramPoint<Real> p)
    {
        if (p.type == DiagramPoint<Real>::DIAG) {
            output << "(" << p.x << ", " << p.y << ", " << 0.5 * (p.x + p.y) << " DIAG)";
        } else {
            output << "(" << p.x << ", " << p.y <<  ")";
        }
        return output;
    }

    template<class Real>
    struct DiagramPointHash {
        size_t operator()(const DiagramPoint<Real>& p) const
        {
            std::size_t seed = 0;
            hash_combine(seed, std::hash<Real>{}(p.x));
            hash_combine(seed, std::hash<Real>{}(p.y));
            hash_combine(seed, std::hash<bool>{}(p.is_diagonal()));
            return seed;
        }
    };

    template<class Real, class Pt>
    struct DistImpl {
        Real operator()(const Pt& a, const Pt& b, const Real p, const int dim)
        {
            Real result = 0.0;
            if (hera::is_infinity(p)) {
                for(int d = 0; d < dim; ++d) {
                    result = std::max(result, std::fabs(a[d] - b[d]));
                }
            } else if (p == 1.0) {
                for(int d = 0; d < dim; ++d) {
                    result += std::fabs(a[d] - b[d]);
                }
            } else {
                assert(p > 1.0);
                for(int d = 0; d < dim; ++d) {
                    result += std::pow(std::fabs(a[d] - b[d]), p);
                }
                result = std::pow(result, 1.0 / p);
            }
            return result;
        }
    };

    template<class Real>
    struct DistImpl<Real, DiagramPoint<Real>> {
        Real operator()(const DiagramPoint<Real>& a, const DiagramPoint<Real>& b, const Real p, const int /*dim */)
        {
            Real result = 0.0;
            if (a.is_diagonal() and b.is_diagonal()) {
                return result;
            } else if (hera::is_infinity(p)) {
                result = std::max(std::fabs(a.getRealX() - b.getRealX()), std::fabs(a.getRealY() - b.getRealY()));
            } else if (p == 1.0) {
                result = std::fabs(a.getRealX() - b.getRealX()) + std::fabs(a.getRealY() - b.getRealY());
            } else {
                assert(p > 1.0);
                result = std::pow(std::pow(std::fabs(a.getRealX() - b.getRealX()), p) + std::pow(std::fabs(a.getRealY() - b.getRealY()), p), 1.0 / p);
            }
            return result;
        }
    };

    template<class R, class Pt>
    inline R dist_lp(const Pt& a, const Pt& b, const R p, const int dim)
    {
        return DistImpl<R, Pt>()(a, b, p, dim);
    }

    template<class Real>
    inline Real dist_l_inf(const DiagramPoint<Real>& a, const DiagramPoint<Real>& b)
    {
        return DistImpl<Real, DiagramPoint<Real>>()(a, b, get_infinity(), 2);
    }
} // namespace hera
#endif
