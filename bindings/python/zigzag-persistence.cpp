#include <vector>
#include <type_traits>

#include <boost/range/adaptor/transformed.hpp>
namespace ba = boost::adaptors;

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <dionysus/row-reduction.h>
#include <dionysus/ordinary-persistence.h>
#include <dionysus/standard-reduction.h>

#include "filtration.h"
#include "zigzag-persistence.h"
#include "diagram.h"

struct Time
{
            Time(float t_, size_t i_, unsigned dim_, bool dir_):
                t(t_), i(i_), dim(dim_), dir(dir_)             {}

    bool    operator<(const Time& other) const
    {
        if (t == other.t)
        {
            if (dir && !other.dir)      // add comes before remove
                return true;
            else if (!dir && other.dir)
                return false;
            else if (dim == other.dim)
                return i < other.i;
            else if (dir)
                return dim < other.dim;
            else // if (!dir)
                return other.dim < dim;
        }
        else
            return t < other.t;
    }

    float       t;
    size_t      i;
    unsigned    dim;
    bool        dir;
};



std::tuple<PyZigzagPersistence, std::vector<PyDiagram>>
zigzag_homology_persistence(const PyFiltration& f, const std::vector<std::vector<float>>& times_, PyZpField::Element prime)
{
    using Index          = PyZigzagPersistence::Index;
    using CellChainEntry = dionysus::ChainEntry<PyZpField, PySimplex>;
    using ChainEntry     = dionysus::ChainEntry<PyZpField, Index>;

    std::vector<Time> times;
    for (size_t i = 0; i < times_.size(); ++i)
    {
        int dim = f[i].dimension();
        bool dir = true;
        for (float t : times_[i])
        {
            times.emplace_back(t, i, dim, dir);
            dir = !dir;
        }
    }
    std::sort(times.begin(), times.end());

    std::vector<PyDiagram> diagrams;
    PyZpField field(prime);
    PyZigzagPersistence persistence(field);
    unsigned op = 0;
    unsigned cell = 0;
    std::vector<unsigned>   cells(f.size(), -1);
    for (auto& tt : times)
    {
        size_t i = tt.i; float t = tt.t; bool dir = tt.dir;

        auto& c = f[i];
        if (dir)
        {
            cells[i] = cell++;

            Index pair = persistence.add(c.boundary(persistence.field()) |
                                                    ba::transformed([&](const CellChainEntry& e)
                                                    {
                                                        auto idx = f.index(e.index());
                                                        return ChainEntry(e.element(), cells[idx]);
                                                    }));

            if (pair != persistence.unpaired())
            {
                auto t_birth = times[pair].t;
                if (t_birth != t)
                {
                    int dim = c.dimension()-1;
                    while (dim+1 > diagrams.size())
                        diagrams.emplace_back();
                    diagrams[dim].emplace_back(t_birth, t, pair);
                }
            }

            ++op;
        } else
        {
            Index pair = persistence.remove(cells[i]);
            cells[i] = -1;
            if (pair != persistence.unpaired())
            {
                auto t_birth = times[pair].t;
                if (t_birth != t)
                {
                    int dim = c.dimension();
                    while (dim+1 > diagrams.size())
                        diagrams.emplace_back();
                    diagrams[dim].emplace_back(t_birth, t, pair);
                }
            }
            ++op;
        }
    }

    // add infinite points
    constexpr float inf = std::numeric_limits<float>::infinity();
    for (auto& birth_idx : persistence.alive_ops())
    {
        auto i_birth   = times[birth_idx].i;
        auto t_birth   = times[birth_idx].t;
        auto dir_birth = times[birth_idx].dir;

        int dim = f[i_birth].dimension();
        if (!dir_birth)     // born on removal
            dim -= 1;
        while (dim+1 > diagrams.size())
            diagrams.emplace_back();
        diagrams[dim].emplace_back(t_birth, inf, birth_idx);
    }

    return std::make_tuple(std::move(persistence), std::move(diagrams));
}

PYBIND11_MAKE_OPAQUE(PyZigzagPersistence::Column);

// This is horribly hacky. In chain.h, we need to output an entry in the
// PyZigzagPersistence::RowMatrix, which is a tuple. So we just provide an
// operator<< to handle it explicitly.
std::ostream& operator<<(std::ostream& out, const PyZigzagPersistence::RowMatrix::Entry::IndexPair& x)
{
    out << std::get<0>(x);
    return out;
}
#include "chain.h"

// FIXME: need to translate chains from times to simplices;
//        fix ZigzagPersistence::add() to take cell_index argument

void init_zigzag_persistence(py::module& m)
{
    using namespace pybind11::literals;
    m.def("zigzag_homology_persistence",   &zigzag_homology_persistence, "filtration"_a, "times"_a, py::arg("prime") = 2,
          R"(
          compute zigzag homology persistence of the filtration with respect to the given times

          Args:
              filtration: an instance of :class:`~dionysus._dionysus.Filtration` with the set of simplices used in the zigzag construction
              times:      a list of lists; the outer list runs parallel with the filtration;
                          the inner list specifies for each simplex when it enters and leaves the zigzag
                          (even entries, starting the indexing from 0, are interpreted as appearance times, odd entires as disappearance)
              prime:      prime modulo which to perform computation

          Returns:
              A pair. The first element is an instance of
              :class:`~dionysus._dionysus.ZigzagPersistence`, which offers access to the cycles
              alive at the end of the zigzag; the second is a list of persistence diagrams.

          )");

    py::class_<PyZigzagPersistence>(m, "ZigzagPersistence", "representation of the current homology basis")
        .def(py::init<PyZpField>())
        .def("__len__",     &PyZigzagPersistence::alive_size,       "number of alive cycles")
        .def("__iter__",    [](const PyZigzagPersistence& zz)
                            {
                                auto cycles = zz.alive_cycles() |
                                                ba::transformed([&zz](PyZigzagPersistence::Index i)
                                                                { return zz.cycle(i); });
                                return py::make_iterator(std::begin(cycles), std::end(cycles));
                            },
                            py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */,
                            "iterator over the alive cycles")
        .def("__repr__",    [](const PyZigzagPersistence& zz)
                            { std::ostringstream oss; oss << "Zigzag persistence with " << zz.alive_size() << " alive cycles"; return oss.str(); })
    ;

    init_chain<PyZigzagPersistence::Column>(m, "ZZ");
}
