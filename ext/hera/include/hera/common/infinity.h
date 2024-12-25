#ifndef HERA_COMMON_INFINITY_H
#define HERA_COMMON_INFINITY_H

namespace hera {

    // to make templates work with Real types that don't have std::numeric_limits<Real>::infinity
    // is used only to denote ell_infinity norm in plane
    template<class Real = double>
    inline constexpr Real get_infinity()
    {
        return Real(-1);
    }

    template<class Real = double>
    inline bool is_infinity(const Real& x)
    {
        return x == get_infinity<Real>();
    }

}

#endif

