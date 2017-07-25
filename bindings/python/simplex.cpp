#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
namespace py = pybind11;

#include "simplex.h"

void init_simplex(py::module& m)
{
    using namespace pybind11::literals;

    py::class_<PySimplex>(m, "Simplex", "an abstract `simplex <https://en.wikipedia.org/wiki/Simplex>`_")
        .def(py::init<>(), "construct empty simplex")
        .def(py::init<std::vector<PySimplex::Vertex>>(), "construct from a list of vertices")
        .def(py::init<std::vector<PySimplex::Vertex>, PySimplex::Data>(), "construct from a list of vertices and data")
        .def("__repr__",        [](const PySimplex& s) { std::ostringstream oss; oss << s << ' ' << s.data(); return oss.str(); })
        .def("boundary",        [](const PySimplex& s) { return py::make_iterator(s.boundary_begin(), s.boundary_end()); },
                                py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */,
                                "returns iterator over the boundary of the simplex")
        .def("dimension",       &PySimplex::dimension, "simplex dimension, one less than cardinality")
        .def("__len__",         &PySimplex::size, "simplex cardinality")
        .def("__getitem__",     &PySimplex::operator[], "access `i`-th vertex", "i"_a)
        .def("__iter__",        [](const PySimplex& s) { return py::make_iterator(s.begin(), s.end()); },
                                py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */,
                                "iterator over the vertices")
        .def("__contains__",    [](const PySimplex& s, PySimplex::Vertex v) { return std::find(s.begin(), s.end(), v) != s.end(); },
                                "test whether the simplex contains given vertex", "v"_a)
        .def("join",            &PySimplex::join, "v"_a, "join a simplex and a vertex")
        .def_property("data",   [](const PySimplex& s) { return s.data(); },
                                [](PySimplex& s, PySimplex::Data x) { s.data() = x; }, "access the data associated to the simplex")
        .def("__hash__",        [](const PySimplex& s) { return hash_value(s); }, "hash simplex")
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self <  py::self)
        .def(py::self >  py::self)
    ;
}
