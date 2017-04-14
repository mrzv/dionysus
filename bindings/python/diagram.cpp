#include <memory>
#include <limits>
#include <iostream>

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "diagram.h"

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
        .def_readwrite("birth",  &Point::birth, "birth value")
        .def_readwrite("death",  &Point::death, "death value")
        .def_readwrite("data",   &Point::data,  "auxiliary data associated to the point (e.g., birth index)")
        .def("__repr__",        [](const Point& p)              { std::ostringstream oss; oss << '(' << p.birth << ',' << p.death << ')'; return oss.str(); })
    ;
}

