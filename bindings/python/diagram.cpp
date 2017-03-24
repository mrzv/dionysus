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
    std::vector<PyDiagram> diagrams;

    for (PyReducedMatrix::Index i = 0; i < m.size(); ++i)
    {
        auto& s = f[i];
        auto  d = s.dimension();

        while (d + 1 > diagrams.size())
            diagrams.emplace_back(PyDiagram());

        auto pair = m.pair(i);
        if (pair == m.unpaired())
        {
            auto  birth = s.data();
            using Value = decltype(birth);
            Value death = std::numeric_limits<Value>::infinity();
            diagrams[d].emplace_back(birth, death, i);
        } else if (pair > i)       // positive
        {
            auto birth = s.data();
            auto death = f[pair].data();

            if (birth != death)         // skip diagonal
                diagrams[d].emplace_back(birth, death, i);
        } // else negative: do nothing
    }

    return diagrams;
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

