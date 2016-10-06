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

    py::class_<PySimplex>(m, "Simplex")
        .def(py::init<>())
        .def(py::init<std::vector<PySimplex::Vertex>>())
        .def(py::init<std::vector<PySimplex::Vertex>, PySimplex::Data>())
        .def("__repr__",        [](const PySimplex& s) { std::ostringstream oss; oss << s << ' ' << s.data(); return oss.str(); })
        .def("boundary",        [](py::object s) { return PySimplexBoundaryIterator(s.cast<const PySimplex&>(), s); })
        .def("dimension",       &PySimplex::dimension)
        .def("__len__",         &PySimplex::size)
        .def("__getitem__",     &PySimplex::operator[])
        .def("__iter__",        [](const PySimplex& s) { return py::make_iterator(s.begin(), s.end()); },
                                py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */)
        .def("__contains__",    [](const PySimplex& s, PySimplex::Vertex v) { return std::find(s.begin(), s.end(), v) != s.end(); })
        .def("join",            &PySimplex::join)
        .def_property("data",   [](const PySimplex& s) { return s.data(); },
                                [](PySimplex& s, PySimplex::Data x) { s.data() = x; })
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self <  py::self)
        .def(py::self >  py::self)
    ;
}
