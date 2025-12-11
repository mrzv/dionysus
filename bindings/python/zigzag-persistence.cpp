#include <limits>
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
        if (t != other.t)
            return t < other.t;

        if (dir != other.dir)
            return dir;

        if (dim != other.dim)
        {
            if (dir)
                return dim < other.dim;
            else
                return other.dim < dim;
        }

        return i < other.i;
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
    std::stable_sort(times.begin(), times.end());

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
                                                        auto idx = f.index(e.index(),i);
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

PyLinkedMultiFiltration
fast_zigzag(const PyFiltration&     f,
            const Times&            times)
{
    int w = -1;
    float inf = std::numeric_limits<float>::infinity();
    PyLinkedMultiFiltration combined;
    combined.push_back(PySimplex({ w }, -inf), 0);

    for (size_t i = 0; i < f.size(); ++i)
    {
        size_t j = 0;
        for (; j < times[i].size(); ++j)
        {
            if (j % 2 == 0)
                combined.push_back(PySimplex(f[i], times[i][j]), combined.size());
            else
                combined.push_back(PySimplex(f[i], times[i][j]).join(w), combined.size() - 1);        // link to the previous appearance
        }

        // if a simplex doesn't get removed, remove it at infinity
        if (j % 2 != 0)
            combined.push_back(PySimplex(f[i], inf).join(w), combined.size() - 1);        // link to the previous appearance
    }

    DataDimCmp base_cmp;
    DataDimCmp cone_cmp(true);
    combined.sort([w,base_cmp,cone_cmp](const PySimplex& x, const PySimplex& y)
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

template<class PyReducedMatrix, class Filtration>
std::vector<std::map<std::string, PyDiagram>>
init_zigzag_diagrams(const PyReducedMatrix& r, const Filtration& f, bool diagonal)
{
    using Index = typename PyReducedMatrix::Index;

    int w = -1;

    std::vector<std::map<std::string, PyDiagram>> result;
    auto result_append = [&result](int dim, std::string type, PySimplex::Data birth, PySimplex::Data death, Index i)
    {
        while (dim >= result.size())
            result.emplace_back();

        result[dim][type].emplace_back(birth,death,i);
    };

    for (Index i = 1; i < r.size(); ++i)
    {
        Index j = r.pair(i);
        if (j < i) continue;        // skip negative

        assert(j != r.unpaired());

        auto i_data = f[i].data();
        auto j_data = f[j].data();

        if (!diagonal && i_data == j_data) continue;

        bool i_cone = f[i].contains(w);
        bool j_cone = f[j].contains(w);

        if (!i_cone && !j_cone)
        {
            // ordinary (closed-open)
            result_append(f[i].dimension(), "co", i_data, j_data, i);
        } else if (i_cone && j_cone)
        {
            // relative (open-closed)
            result_append(f[i].dimension() - 1, "oc", j_data, i_data, i);
        } else
        {
            assert(!i_cone && j_cone);
            if (i_data > j_data)        // TODO: can we check this non-numerically
            {
                // extended (open-open)
                result_append(f[i].dimension() - 1, "oo", j_data, i_data, i);
            } else
            {
                // extended (closed-closed)
                result_append(f[i].dimension(), "cc", i_data, j_data, i);
            }
        }
    }

    return result;
}

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

    m.def("fast_zigzag",    &fast_zigzag, "filtration"_a, "times"_a,
          "Build the cone to compute extended persistence equivalent to the given zigzag.");
    m.def("init_zigzag_diagrams",    &init_zigzag_diagrams<PyReducedMatrix,PyLinkedMultiFiltration>,
          "r"_a, "f"_a, "diagonal"_a = false,
          "Given the cone `f` and its reduced matrix `r`, initialize zigzag diagrams.");
    m.def("init_zigzag_diagrams",    &init_zigzag_diagrams<PyReducedMatrixWithV,PyLinkedMultiFiltration>,
          "r"_a, "f"_a, "diagonal"_a = false,
          "Given the cone `f` and its reduced matrix `r`, initialize zigzag diagrams.");
    m.def("init_zigzag_diagrams",    &init_zigzag_diagrams<PyReducedMatrixNoNegative,PyLinkedMultiFiltration>,
          "r"_a, "f"_a, "diagonal"_a = false,
          "Given the cone `f` and its reduced matrix `r`, initialize zigzag diagrams.");
    m.def("init_zigzag_diagrams",    &init_zigzag_diagrams<PyReducedMatrixNoNegativeWithV,PyLinkedMultiFiltration>,
          "r"_a, "f"_a, "diagonal"_a = false,
          "Given the cone `f` and its reduced matrix `r`, initialize zigzag diagrams.");
}
