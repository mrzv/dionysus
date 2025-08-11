#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
namespace py = pybind11;

#include "simplex.h"
#include "filtration.h"

template<class PyFiltration>
void export_filtration(py::class_<PyFiltration>& cls)
{
    using namespace pybind11::literals;
    cls
        .def(py::init<>(),      "initialize empty filtration")
        .def("__len__",         &PyFiltration::size,        "size of the filtration")
        .def("__getitem__",     &PyFiltration::operator[],  "access the simplex at the given index")
        .def("__setitem__",     &PyFiltration::replace,     "replace the simplex at the given index")
        .def("index",           [](PyFiltration* f, const PySimplex& s) { return f->index(s,f->size()); }, "s"_a, "find an ordered index of a simplex s in the filtration")
        .def("index",           &PyFiltration::index,       "s"_a, "i"_a, "find the ordered index of a simplex s (no later than i) in the filtration")
        .def("__contains__",    &PyFiltration::contains,    "test whether filtration contains the simplex")
        .def("__iter__",        [](const PyFiltration& f) { return py::make_iterator(f.begin(), f.end()); },
                                py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */,
                                "iterate over the simplices in sorted order")
        .def("sort",            [](PyFiltration& f, bool reverse) { f.sort(DataDimCmp(reverse)); },
                                "reverse"_a = false,
                                "sort the filtration with respect to data, breaking ties using dimension, and then lexicographically")
        .def("sort",            [](PyFiltration& f, std::function<int (const PySimplex&, const PySimplex&)> cmp, bool reverse)
                                {
                                    std::cerr << "Warning: cmp-based sort is deprecated (to match Python 3); use key-based sort instead" << std::endl;
                                    f.sort([cmp,reverse](const PySimplex& s1, const PySimplex& s2) { return reverse ? cmp(s1, s2) > 0 : cmp(s1, s2) < 0; });
                                },
                                "cmp"_a, "reverse"_a = false,
                                "sort the filtration with respect to the given functor")
        .def("sort",            [](PyFiltration& f, std::function<py::tuple (const PySimplex& s)> key, bool reverse)
                                {
                                    f.sort([key,reverse](const PySimplex& s1, const PySimplex& s2)
                                           {
                                               return reverse ? key(s1) > key(s2) : key(s1) < key(s2);
                                           });
                                },
                                py::kw_only(), "key"_a, "reverse"_a = false,
                                "sort the filtration with respect to the given key")
        .def("rearrange",       &PyFiltration::rearrange, "indices"_a, "rearrange simplices into the given order")
        .def("__repr__",        [](const PyFiltration& f) { std::ostringstream oss; oss << "Filtration with " << f.size() << " simplices"; return oss.str(); })
    ;
}

template<class PyFiltration>
void export_ordinary(py::module& m, std::string name)
{
    using namespace pybind11::literals;
    py::class_<PyFiltration> cls(m, name.c_str(), "store an ordered sequence of simplices, providing lookup");
    cls
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
    ;
    export_filtration(cls);
}

template<class PyFiltration>
void export_linked(py::module& m, std::string name)
{
    using namespace pybind11::literals;
    py::class_<PyFiltration> cls(m, name.c_str(), "store an ordered sequence of simplices, providing lookup");
    cls
        .def(py::init([](const std::vector<std::tuple<PySimplex, size_t>>& simplices)
                      {
                          PyFiltration* f = new PyFiltration;
                          for (auto& sl : simplices)
                              f->emplace_back(std::get<1>(sl), std::get<0>(sl));
                          return f;
                      }), "initialize filtration from a list of simplices")
        .def(py::init([](const std::vector<std::tuple<std::vector<PySimplex::Vertex>, size_t>>& simplices)
                      {
                          PyFiltration* f = new PyFiltration;
                          for (auto& sl : simplices)
                              f->emplace_back(std::get<1>(sl), std::get<0>(sl));
                          return f;
                      }), "initialize filtration from a list of lists")
        .def(py::init([](const std::vector<std::tuple<std::vector<PySimplex::Vertex>, PySimplex::Data, size_t>>& simplices)
                      {
                          PyFiltration* f = new PyFiltration;
                          for (auto& sl : simplices)
                              f->emplace_back(std::get<2>(sl), std::get<0>(sl), std::get<1>(sl));
                          return f;
                      }), "initialize filtration from a list of tuples of vertices and values")
        .def("append",  [](PyFiltration* f, const PySimplex& s, size_t l) { f->push_back(s,l); },  "s"_a, "l"_a, "append simplex to the filtration")
    ;
    export_filtration(cls);
}

void init_filtration(py::module& m)
{
    export_ordinary<PyFiltration>(m, "Filtration");
    export_ordinary<PyMultiFiltration>(m, "MultiFiltration");
    export_linked<PyLinkedMultiFiltration>(m, "LinkedMultiFiltration");
}
