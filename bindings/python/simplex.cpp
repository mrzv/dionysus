#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
namespace py = pybind11;

#include "simplex.h"

// Need special handling of the boundary iterator
struct PySimplexBoundaryIterator
{
    PySimplexBoundaryIterator(const PySimplex& s, py::object ref):
        s(s), ref(ref), it(s.boundary_begin())  {}

    PySimplex next()
    {
        if (it == s.boundary_end())
            throw py::stop_iteration();
        return *it++;
    }

    const PySimplex& s;
    py::object ref;
    PySimplex::BoundaryIterator it;
};

void init_simplex(py::module& m)
{
    py::class_<PySimplexBoundaryIterator>(m, "SimplexBoundaryIterator")
        .def("__iter__", [](PySimplexBoundaryIterator& it) -> PySimplexBoundaryIterator& { return it; })
        .def("__next__", &PySimplexBoundaryIterator::next);

    py::class_<PySimplex>(m, "Simplex", "an abstract `simplex <https://en.wikipedia.org/wiki/Simplex>`_")
        .def(py::init<>(), "construct empty simplex")
        .def(py::init<std::vector<PySimplex::Vertex>>(), "construct from a list of vertices")
        .def(py::init<std::vector<PySimplex::Vertex>, PySimplex::Data>(), "construct from a list of vertices and data")
        .def("__repr__",        [](const PySimplex& s) { std::ostringstream oss; oss << s << ' ' << s.data(); return oss.str(); })
        .def("boundary",        [](py::object s) { return PySimplexBoundaryIterator(s.cast<const PySimplex&>(), s); },
                                "returns iterator over the boundary of the simplex")
        .def("dimension",       &PySimplex::dimension, "simplex dimension, one less than cardinality")
        .def("__len__",         &PySimplex::size, "simplex cardinality")
        .def("__getitem__",     &PySimplex::operator[], "access $i$-th vertex")
        .def("__iter__",        [](const PySimplex& s) { return py::make_iterator(s.begin(), s.end()); },
                                py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */,
                                "iterator over the vertices")
        .def("__contains__",    [](const PySimplex& s, PySimplex::Vertex v) { return std::find(s.begin(), s.end(), v) != s.end(); },
                                "test whether the simplex contains given vertex")
        .def("join",            &PySimplex::join, "join two simplices, i.e., take the union of their vertices")
        .def_property("data",   [](const PySimplex& s) { return s.data(); },
                                [](PySimplex& s, PySimplex::Data x) { s.data() = x; }, "access the data associated to the simplex")
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self <  py::self)
        .def(py::self >  py::self)
    ;
}
