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
        .def(py::init([](const std::vector<PySimplex>& simplices)
                      {
                          PyFiltration* f = new PyFiltration;
                          for (auto& s : simplices)
                              f->emplace_back(s);
                          return f;
                      }), "initialize filtration from a list of simplices")
        .def(py::init([](const std::vector<std::vector<PySimplex::Vertex>>& simplices)
                      {
                          PyFiltration* f = new PyFiltration;
                          for (auto& s : simplices)
                              f->emplace_back(s);
                          return f;
                      }), "initialize filtration from a list of lists")
        .def(py::init([](const std::vector<std::tuple<std::vector<PySimplex::Vertex>, PySimplex::Data>>& simplices)
                      {
                          PyFiltration* f = new PyFiltration;
                          for (auto& s : simplices)
                              f->emplace_back(std::get<0>(s), std::get<1>(s));
                          return f;
                      }), "initialize filtration from a list of tuples of vertices and values")
        .def("append",          [](PyFiltration* f, const PySimplex& s) { f->push_back(s); },  "s"_a, "append simplex to the filtration")
        .def("add",             [](PyFiltration* f, const PySimplex& s) { return f->add(s); }, "s"_a,
                                "append simplex to the filtration, if not already in the filtration; either way return the index of the simplex")
        .def("__len__",         &PyFiltration::size,        "size of the filtration")
        .def("__getitem__",     &PyFiltration::operator[],  "access the simplex at the given index")
        .def("__setitem__",     &PyFiltration::replace,     "replace the simplex at the given index")
        .def("index",           &PyFiltration::index,       "s"_a, "find the ordered index of a simplex in the filtration")
        .def("__contains__",    &PyFiltration::contains,    "test whether filtration contains the simplex")
        .def("__iter__",        [](const PyFiltration& f) { return py::make_iterator(f.begin(), f.end()); },
                                py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */,
                                "iterate over the simplices in sorted order")
        .def("sort",            [](PyFiltration& f, bool reverse) { f.sort(DataDimCmp(reverse)); },
                                "reverse"_a = false,
                                "sort the filtration with respect to data, breaking ties using dimension, and then lexicographically")
        .def("sort",            [](PyFiltration& f, std::function<int (const PySimplex&, const PySimplex&)> cmp, bool reverse)
                                { f.sort([cmp,reverse](const PySimplex& s1, const PySimplex& s2) { return reverse ? cmp(s1, s2) > 0 : cmp(s1, s2) < 0; }); },
                                "cmp"_a, "reverse"_a = false,
                                "sort the filtration with respect to the given functor")
        .def("rearrange",       &PyFiltration::rearrange, "indices"_a, "rearrange simplices into the given order")
        .def("__repr__",        [](const PyFiltration& f) { std::ostringstream oss; oss << "Filtration with " << f.size() << " simplices"; return oss.str(); })
    ;
}
