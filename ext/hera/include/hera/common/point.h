#ifndef HERA_POINT_H
#define HERA_POINT_H

#include <ostream>

namespace hera {

    template<class Real = double>
    struct Point {
        Real x, y;

        bool operator==(const Point<Real>& other) const { return ((this->x == other.x) and (this->y == other.y)); }
        bool operator!=(const Point<Real>& other) const { return !(*this == other); }

        Point(Real ax, Real ay) :x(ax), y(ay) { }
        Point() :x(0.0), y(0.0) { }

        template<class R>
        friend std::ostream& operator<<(std::ostream& output, const Point<R>& p)
        {
            output << "(" << p.x << ", " << p.y << ")";
            return output;
        }
    };

    template<class Real>
    inline Real sqr_dist(const Point<Real>& a, const Point<Real>& b)
    {
        return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
    }

    template<class Real>
    inline Real dist(const Point<Real>& a, const Point<Real>& b)
    {
        return sqrt(sqr_dist(a, b));
    }

} // namespace hera

#endif
