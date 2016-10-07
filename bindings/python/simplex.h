#pragma once

#include <dionysus/simplex.h>

using PySimplex = dionysus::Simplex<int, float>;

inline bool data_dim_cmp(const PySimplex& x, const PySimplex& y)
{
    return x.data() < y.data() || (x.data() == y.data() && x < y);      // x < y compares dimension first and then compares lexicographically
}
