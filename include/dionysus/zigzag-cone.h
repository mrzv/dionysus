#ifndef DIONYSUS_ZIGZAG_CONE_H
#define DIONYSUS_ZIGZAG_CONE_H

#include <cassert>
#include <limits>
#include <map>
#include <string>
#include <vector>

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

template<class ReducedMatrix, class Filtration, class Diagram>
std::vector<std::map<std::string, Diagram>> init_zigzag_diagrams(const ReducedMatrix& r,
                                                                 const Filtration& f,
                                                                 bool diagonal,
                                                                 typename Filtration::Cell::Vertex cone_vertex = -1)
{
    using Index = typename ReducedMatrix::Index;
    using Data = typename Filtration::Cell::Data;

    auto w = cone_vertex;
    std::vector<std::map<std::string, Diagram>> result;
    auto result_append = [&result](int dim, std::string type, Data birth, Data death, Index i)
    {
        while (dim >= result.size())
            result.emplace_back();

        result[dim][type].emplace_back(birth, death, i);
    };

    for (Index i = 1; i < r.size(); ++i)
    {
        Index j = r.pair(i);
        if (j < i) continue;

        assert(j != r.unpaired());

        auto i_data = f[i].data();
        auto j_data = f[j].data();

        if (!diagonal && i_data == j_data) continue;

        bool i_cone = f[i].contains(w);
        bool j_cone = f[j].contains(w);

        if (!i_cone && !j_cone)
        {
            result_append(f[i].dimension(), "co", i_data, j_data, i);
        } else if (i_cone && j_cone)
        {
            result_append(f[i].dimension() - 1, "oc", j_data, i_data, i);
        } else
        {
            assert(!i_cone && j_cone);
            if (i_data > j_data)
                result_append(f[i].dimension() - 1, "oo", j_data, i_data, i);
            else
                result_append(f[i].dimension(), "cc", i_data, j_data, i);
        }
    }

    return result;
}

}

#endif
