#include <vector>
#include <type_traits>

#include <boost/range/adaptor/transformed.hpp>
namespace ba = boost::adaptors;

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
namespace py = pybind11;

#include <dionysus/row-reduction.h>
#include <dionysus/ordinary-persistence.h>
#include <dionysus/standard-reduction.h>

#include "filtration.h"
#include "persistence.h"                // to get access to PyReducedMatrix::Chain
#include "zigzag-persistence.h"
#include "diagram.h"
#include "progress.h"

PYBIND11_MAKE_OPAQUE(PyReducedMatrix::Chain);      // we want to provide our own binding for Chain

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

class PyTimeIndexMap
{
    public:
        using ZZIndex           = PyZigzagPersistence::Index;
        using FIndex            = size_t;
        using Map               = std::unordered_map<ZZIndex, FIndex>;
        using const_iterator    = Map::const_iterator;

        void                    set(const ZZIndex& x, const FIndex& y)  { m_[x] = y; }
        void                    remove(const ZZIndex& x)                { m_.erase(x); }

        FIndex                  operator[](const ZZIndex& x) const      { return m_.find(x)->second; }
        size_t                  size() const                            { return m_.size(); }

        const_iterator          begin() const                           { return m_.begin(); }
        const_iterator          end() const                             { return m_.end(); }

    private:
        Map     m_;
};

using Times     = std::vector<std::vector<float>>;
using Callback  = std::function<void(size_t, float, bool, const PyZigzagPersistence*, const PyTimeIndexMap*)>;

std::tuple<PyZigzagPersistence, std::vector<PyDiagram>, PyTimeIndexMap>
zigzag_homology_persistence(const PyFiltration&     f,
                            const Times&            times_,
                            PyZpField::Element      prime,
                            const Callback&         callback,
                            bool                    show_progress)
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

    std::unique_ptr<Progress> progress(new NoProgress);
    if (show_progress)
        progress = std::unique_ptr<Progress>(new ShowProgress(times.size()));

    std::vector<PyDiagram> diagrams;
    PyZpField field(prime);
    PyZigzagPersistence persistence(field);
    unsigned op = 0;
    unsigned cell = 0;
    std::vector<unsigned>   cells(f.size(), -1);
    PyTimeIndexMap          cells_inv_;
    for (auto& tt : times)
    {
        (*progress)();

        size_t i = tt.i; float t = tt.t; bool dir = tt.dir;

        auto& c = f[i];
        if (dir)
        {
            cells_inv_.set(cell, i);
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
            cells_inv_.remove(cells[i]);
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

        callback(i,t,dir,&persistence,&cells_inv_);
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

    return std::make_tuple(std::move(persistence), std::move(diagrams), std::move(cells_inv_));
}

// custom iterator that preserves right filtration order
struct PyZZAliveCycleIterator
{
            PyZZAliveCycleIterator(const PyZigzagPersistence& zz_, py::object ref_):
                zz(zz_), ref(ref_)
    {
        for (PyZigzagPersistence::Index x : zz.alive_cycles())
            indices.push_back(x);
        std::sort(indices.begin(), indices.end());
    }

    PyReducedMatrix::Chain  next()
    {
        if (idx == indices.size())
            throw py::stop_iteration();

        PyReducedMatrix::Chain c;
        for (auto& x : zz.cycle(indices[idx]))
            c.emplace_back(x.element(), std::get<0>(x.index()));
        ++idx;
        return c;
    }

    const PyZigzagPersistence&                  zz;
    py::object                                  ref;
    size_t                                      idx = 0;
    std::vector<PyZigzagPersistence::Index>     indices;
};

#include "chain.h"

void init_zigzag_persistence(py::module& m)
{
    using namespace pybind11::literals;
    m.def("zigzag_homology_persistence",   &zigzag_homology_persistence, "filtration"_a, "times"_a, "prime"_a = 2,
                                                                         "callback"_a = Callback([](size_t, float, bool, const PyZigzagPersistence*, const PyTimeIndexMap*){}),
                                                                         "progress"_a = false,
          R"(
          compute zigzag homology persistence of the filtration with respect to the given times

          Args:
              filtration: an instance of :class:`~dionysus._dionysus.Filtration` with the set of simplices used in the zigzag construction
              times:      a list of lists; the outer list runs parallel with the filtration;
                          the inner list specifies for each simplex when it enters and leaves the zigzag
                          (even entries, starting the indexing from 0, are interpreted as appearance times, odd entires as disappearance)
              prime:      prime modulo which to perform computation
              callback:   function to call after every step in the zigzag; it gets arguments `(i,t,d,zz,cells)`,
                          where `i` is the index of the simplex being added or removed, `t` is the time,
                          `d` is the "direction" (`True` if the simplex is being added, `False` if it`s being removed),
                          `zz` is the current state of the :class:`~dionysus._dionysus.ZigzagPersistence`,
                          `cells` is the map from the internal indices of the zigzag representation to the filtration indices,
              progress:   show a progress bar.

          Returns:
              A triple. The first element is an instance of
              :class:`~dionysus._dionysus.ZigzagPersistence`, which offers access to the cycles
              alive at the end of the zigzag; the second is a list of
              persistence diagrams; the third is an instance of
              :class:`~dionysus._dionysus.TimeIndexMap` for translating cycles
              from the internal representation to filtration indices.
          )");

   py::class_<PyZZAliveCycleIterator>(m, "ZZAliveCycleIterator")
        .def("__iter__", [](PyZZAliveCycleIterator& it) -> PyZZAliveCycleIterator& { return it; })
        .def("__next__", &PyZZAliveCycleIterator::next);

    py::class_<PyZigzagPersistence>(m, "ZigzagPersistence", "representation of the current homology basis")
        .def(py::init<PyZpField>())
        .def("__len__",     &PyZigzagPersistence::alive_size,       "number of alive cycles")
        .def("__iter__",    [](py::object zz)   { return PyZZAliveCycleIterator(zz.cast<const PyZigzagPersistence&>(), zz); }, "iterator over the alive cycles")
        .def("__repr__",    [](const PyZigzagPersistence& zz)
                            { std::ostringstream oss; oss << "Zigzag persistence with " << zz.alive_size() << " alive cycles"; return oss.str(); })
    ;

    py::class_<PyTimeIndexMap>(m, "TimeIndexMap", "map from the internal representation of zigzag persistence to filtration indices")
        .def("__len__",     &PyTimeIndexMap::size,                      "size of the map")
        .def("__getitem__", [](const PyTimeIndexMap& m, size_t i) { return m[i]; }, "access the filtration index of the given internal index")
        .def("__iter__",    [](const PyTimeIndexMap& m) { return py::make_iterator(m.begin(), m.end()); },
                                py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */,
                                "iterate over the entries of the map")
    ;
}
