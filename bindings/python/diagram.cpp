#include <memory>
#include <limits>
#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "diagram.h"
#include "persistence.h"
#include "filtration.h"

std::vector<PyDiagram>
py_init_diagrams(const PyReducedMatrix& m, const PyFiltration& f)
{
    return init_diagrams(m, f,
                         [](const PySimplex& s)                     { return s.data(); },       // value
                         [](PyReducedMatrix::Index i) -> size_t     { return i; });             // data
}

void init_diagram(py::module& m)
{
    using namespace pybind11::literals;
    py::class_<PyDiagram>(m, "Diagram", "persistence diagram")
        .def(py::init<>(),      "initialize empty filtration")
        .def("append",          &PyDiagram::push_back, "p"_a,   "append point to the diagram")
        .def("__len__",         &PyDiagram::size,               "size of the diagram")
        .def("__iter__",        [](const PyDiagram& dgm) { return py::make_iterator(dgm.begin(), dgm.end()); },
                                py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */,
                                "iterate over the points of the diagram")
        .def("__repr__",        [](const PyDiagram& dgm)        { std::ostringstream oss; oss << "Diagram with " << dgm.size() << " points"; return oss.str(); })
    ;

    using Point = PyDiagram::Point;
    py::class_<Point>(m, "DiagramPoint", "persistence diagram point")
        .def_readwrite("birth",  &Point::birth)
        .def_readwrite("death",  &Point::death)
        .def_readwrite("data",   &Point::data)
        .def("__repr__",        [](const Point& p)              { std::ostringstream oss; oss << '(' << p.birth << ',' << p.death << ')'; return oss.str(); })
    ;

    m.def("init_diagrams",      &py_init_diagrams,  "m"_a, "f"_a,  "initialize diagrams from reduced matrix and filtration");
}

