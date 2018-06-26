#pragma once

#include <dionysus/simplex.h>

using PySimplex = dionysus::Simplex<int, float>;

struct DataDimCmp
{
            DataDimCmp(bool reverse_ = false):
                reverse(reverse_)       {}

    bool    operator()(const PySimplex& x, const PySimplex& y) const
    {
        bool res = reverse ? x.data() > y.data() : x.data() < y.data();
        return res || (x.data() == y.data() && x < y);      // x < y compares dimension first and then compares lexicographically
    }

    bool reverse;
};
