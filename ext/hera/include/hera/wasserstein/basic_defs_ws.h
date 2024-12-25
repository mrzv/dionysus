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

#ifndef BASIC_DEFS_WS_H
#define BASIC_DEFS_WS_H

#include <vector>
#include <math.h>
#include <cstddef>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <iomanip>
#include <locale>
#include <cassert>
#include <limits>
#include <ostream>
#include <typeinfo>

#ifdef _WIN32
#include <ciso646>
#endif

#include <hera/common.h>
#include <hera/dnn/geometry/euclidean-dynamic.h>
#include "def_debug_ws.h"
#include "auction_params.h"

namespace hera
{

template<class Real = double>
inline bool is_p_valid_norm(const Real& p)
{
    return is_infinity<Real>(p) or p >= Real(1);
}

namespace ws
{

    using IdxType = int;

    constexpr size_t k_invalid_index = std::numeric_limits<IdxType>::max();

    template<class Real = double>
    using IdxValPair = std::pair<IdxType, Real>;

    template<class R>
    inline std::ostream& operator<<(std::ostream& output, const IdxValPair<R> p)
    {
        output << "(" << p.first << ", " << p.second << ")";
        return output;
    }


    // TODO
    template<class Real, typename DiagPointContainer>
    inline double getFurthestDistance3Approx(DiagPointContainer& A, DiagPointContainer& B, const Real p)
    {
        int dim = 2;
        Real result { 0.0 };
        DiagramPoint<Real> begA = *(A.begin());
        DiagramPoint<Real> optB = *(B.begin());
        for(const auto& pointB : B) {
            if (dist_lp(begA, pointB, p, dim) > result) {
                result = dist_lp(begA, pointB, p, dim);
                optB = pointB;
            }
        }
        for(const auto& pointA : A) {
            if (dist_lp(pointA, optB, p, dim) > result) {
                result = dist_lp(pointA, optB, p, dim);
            }
        }
        return result;
    }

    template<class Real>
    inline Real getFurthestDistance3Approx_pg(const hera::ws::dnn::DynamicPointVector<Real>& A, const hera::ws::dnn::DynamicPointVector<Real>& B, const Real p, const int dim)
    {
        Real result { 0.0 };
        int opt_b_idx = 0;
        for(size_t b_idx = 0; b_idx < B.size(); ++b_idx) {
            if (dist_lp(A[0], B[b_idx], p, dim) > result) {
                result = dist_lp(A[0], B[b_idx], p, dim);
                opt_b_idx = b_idx;
            }
        }

        for(size_t a_idx = 0; a_idx < A.size(); ++a_idx) {
            result = std::max(result,  dist_lp(A[a_idx], B[opt_b_idx], p, dim));
        }

        return result;
    }

    template<class Container>
    inline std::string format_container_to_log(const Container& cont);

    template<class Real, class IndexContainer>
    inline std::string format_point_set_to_log(const IndexContainer& indices, const std::vector<DiagramPoint<Real>>& points);

    template<class T>
    inline std::string format_int(T i);

} // ws
} // hera



#include "basic_defs_ws.hpp"


#endif
