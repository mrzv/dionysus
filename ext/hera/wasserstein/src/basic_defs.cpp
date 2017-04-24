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

#include <algorithm>
#include <cfloat>
#include <set>
#include "basic_defs_ws.h"

#ifndef FOR_R_TDA
#include <iostream>
#endif

namespace geom_ws {
// Point

bool Point::operator==(const Point& other) const
{
    return ((this->x == other.x) and (this->y == other.y));
}

bool Point::operator!=(const Point& other) const
{
    return !(*this == other);
}


#ifndef FOR_R_TDA
std::ostream& operator<<(std::ostream& output, const Point p)
{
    output << "(" << p.x << ", " << p.y << ")";
    return output;
}
#endif

double sqrDist(const Point& a, const Point& b)
{
    return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}

double dist(const Point& a, const Point& b)
{
    return sqrt(sqrDist(a, b));
}

// DiagramPoint

// compute l-inf distance between two diagram points
double distLInf(const DiagramPoint& a, const DiagramPoint& b)
{
    if (a.isDiagonal() and b.isDiagonal()) {
        return 0.0;
    } else {
        return std::max(fabs(a.getRealX() - b.getRealX()), fabs(a.getRealY() - b.getRealY()));
    }
}

double distLp(const DiagramPoint& a, const DiagramPoint& b, const double p)
{
    // infinity: special case
    if ( std::isinf(p) )
        return distLInf(a, b);

    // check p
    assert( p >= 1.0 );

    // avoid calling pow function
    if ( p == 1.0 ) {
        if ( a.isNormal() or b.isNormal() ) {
            // distance between normal points is a usual l-inf distance
            return fabs(a.getRealX() - b.getRealX()) + fabs(a.getRealY() - b.getRealY());
        } else 
            return 0.0;
    } 

    if ( a.isNormal() or b.isNormal() ) {
        // distance between normal points is a usual l-inf distance
        return std::pow(std::pow(fabs(a.getRealX() - b.getRealX()), p) + std::pow(fabs(a.getRealY() - b.getRealY()), p), 1.0/p );
    } else 
        return 0.0;
}


double DiagramPoint::persistenceLp(const double p) const
{
    if (isDiagonal())
        return 0.0;
    else {
        double u { 0.5 * (getRealY() + getRealX()) };
        DiagramPoint a_proj(u, u, DiagramPoint::DIAG);
        return distLp(*this, a_proj, p);
    }


}


#ifndef FOR_R_TDA
std::ostream& operator<<(std::ostream& output, const DiagramPoint p)
{
    if ( p.type == DiagramPoint::DIAG ) {
        output << "(" << p.x << ", " << p.y << ", " <<  0.5 * (p.x + p.y) << " DIAG )";
    } else {
        output << "(" << p.x << ", " << p.y << ", " << " NORMAL)";
    }
    return output;
}
#endif

DiagramPoint::DiagramPoint(double xx, double yy, Type ttype) : 
    x(xx),
    y(yy),
    type(ttype)
{
    //if ( yy < xx )
        //throw "Point is below the diagonal";
    //if ( yy == xx and ttype != DiagramPoint::DIAG)
        //throw "Point on the main diagonal must have DIAG type";
}

double DiagramPoint::getRealX() const
{
    if (isNormal())
        return x;
    else 
        return 0.5 * ( x + y);
}

double DiagramPoint::getRealY() const
{
    if (isNormal())
        return y;
    else 
        return 0.5 * ( x + y);
}

} // end of namespace geom_ws
