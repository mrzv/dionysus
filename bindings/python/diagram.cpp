#include <memory>
#include <limits>
#include <iostream>
#include <cmath>

#include <boost/functional/hash.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
namespace py = pybind11;

#include "diagram.h"

void init_diagram(py::module& m)
{
    using namespace pybind11::literals;
    py::class_<PyDiagram>(m, "Diagram", "persistence diagram")
        .def(py::init<>(),      "initialize empty diagram")
        .def(py::init([](const std::vector<std::tuple<PySimplex::Data, PySimplex::Data>>& pts)
                      {
                          PyDiagram* dgm = new PyDiagram;
                          for (auto& pt : pts)
                              dgm->emplace_back(std::get<0>(pt), std::get<1>(pt), 0);
                          return dgm;
                      }), "initialize diagram from a list of (birth,death) points")
        .def("append",          &PyDiagram::push_back, "p"_a,   "append point to the diagram")
        .def("__getitem__",     &PyDiagram::operator[],         "access `i`-th point")
        .def("__len__",         &PyDiagram::size,               "size of the diagram")
        .def("__iter__",        [](const PyDiagram& dgm) { return py::make_iterator(dgm.begin(), dgm.end()); },
                                py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */,
                                "iterate over the points of the diagram")
        .def("__repr__",        [](const PyDiagram& dgm)        { std::ostringstream oss; oss << "Diagram with " << dgm.size() << " points"; return oss.str(); })
    ;

    using Point = PyDiagram::Point;
    py::class_<Point>(m, "DiagramPoint", "persistence diagram point")
        .def_property_readonly("birth",  &Point::birth, "birth value")
        .def_property_readonly("death",  &Point::death, "death value")
        .def_readwrite("data",           &Point::data,  "auxiliary data associated to the point (e.g., birth index)")
        .def(py::self == py::self)
        .def("__hash__",        [](const Point& p)              { return boost::hash<Point::Parent>()(p); },        "hash of the point")
        .def("__repr__",        [](const Point& p)              { std::ostringstream oss; oss << '(' << p.birth() << ',' << p.death() << ')'; return oss.str(); })
    ;

    m.def("log", [](const PyDiagram& dgm)
                 {
                    PyDiagram result;
                    for (auto& pt : dgm)
                        result.emplace_back(std::log(pt.birth()), std::log(pt.death()), pt.data);
                    return result;
                 }, "dgm"_a, "take log of persistence diagram");
}

