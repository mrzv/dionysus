#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
namespace py = pybind11;

#include "simplex.h"
#include "filtration.h"

void init_filtration(py::module& m)
{
    using namespace pybind11::literals;
    py::class_<PyFiltration>(m, "Filtration", "store an ordered sequence of simplices, providing lookup")
        .def(py::init<>(),      "initialize empty filtration")
        .def(py::init<std::vector<PySimplex>>(),      "initialize filtration from a list")
        .def("append",          [](PyFiltration* f, const PySimplex& s) { f->push_back(s); },  "s"_a, "append simplex to the filtration")
        .def("add",             [](PyFiltration* f, const PySimplex& s) { return f->add(s); }, "s"_a,
                                "append simplex to the filtration, if not already in the filtration; either way return the index of the simplex")
        .def("__len__",         &PyFiltration::size,        "size of the filtration")
        .def("__getitem__",     &PyFiltration::operator[],  "access the simplex at the given index")
        .def("index",           &PyFiltration::index,       "s"_a, "find the ordered index of a simplex in the filtration")
        .def("__contains__",    &PyFiltration::contains,    "test whether filtration contains the simplex")
        .def("__iter__",        [](const PyFiltration& f) { return py::make_iterator(f.begin(), f.end()); },
                                py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */,
                                "iterate over the simplices in sorted order")
        .def("sort",            [](PyFiltration& f) { f.sort(data_dim_cmp); },
                                "sort the filtration with respect to data, breaking ties using dimension, and then lexicographically")
        .def("sort",            [](PyFiltration& f, std::function<int (const PySimplex&, const PySimplex&)> cmp)
                                { f.sort([cmp](const PySimplex& s1, const PySimplex& s2) { return cmp(s1, s2) < 0; }); },
                                "cmp"_a,
                                "sort the filtration with respect to the given functor")
        .def("__repr__",        [](const PyFiltration& f) { std::ostringstream oss; oss << "Filtration with " << f.size() << " simplices"; return oss.str(); })
    ;
}
