#ifndef DIONYSUS_ZIGZAG_CONE_H
#define DIONYSUS_ZIGZAG_CONE_H

#include <limits>

namespace dionysus
{

template<class LinkedFiltration, class Filtration, class Times, class BaseCompare, class ConeCompare>
LinkedFiltration make_zigzag_cone(const Filtration& f,
                                  const Times& times,
                                  BaseCompare base_cmp,
                                  ConeCompare cone_cmp,
                                  typename LinkedFiltration::Cell::Vertex cone_vertex = -1)
{
    using Cell = typename LinkedFiltration::Cell;
    using Data = typename Cell::Data;

    auto w = cone_vertex;
    auto inf = std::numeric_limits<Data>::infinity();
    LinkedFiltration combined;
    combined.push_back(Cell({ w }, -inf), 0);

    for (size_t i = 0; i < f.size(); ++i)
    {
        size_t j = 0;
        for (; j < times[i].size(); ++j)
        {
            if (j % 2 == 0)
                combined.push_back(Cell(f[i], times[i][j]), combined.size());
            else
                combined.push_back(Cell(f[i], times[i][j]).join(w), combined.size() - 1);
        }

        // If a simplex does not get removed, remove it at infinity.
        if (j % 2 != 0)
            combined.push_back(Cell(f[i], inf).join(w), combined.size() - 1);
    }

    combined.sort([w,base_cmp,cone_cmp](const Cell& x, const Cell& y)
                  {
                      bool x_cone = x.contains(w);
                      bool y_cone = y.contains(w);

                      if (x_cone && x.dimension() == 0)
                          return true;
                      if (y_cone && y.dimension() == 0)
                          return false;

                      if (!x_cone && y_cone) return true;
                      if (x_cone && !y_cone) return false;

                      if (!x_cone)
                          return base_cmp(x,y);
                      else
                          return cone_cmp(x,y);
                  });

    return combined;
}

}

#endif
