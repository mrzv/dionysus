#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "simplex.h"
#include "filtration.h"

void init_filtration(py::module& m)
{
    auto data_dim_cmp = [](const PySimplex& x, const PySimplex& y)
                        {
                            return x.data() < y.data() || (x.data() == y.data() && x < y);      // x < y compares dimension first and then compares lexicographically
                        };

    py::class_<PyFiltration>(m, "Filtration")
        .def(py::init<>())
        .def("add",             [](PyFiltration* f, const PySimplex& s) { f->push_back(s); })
        .def("__len__",         &PyFiltration::size)
        .def("__getitem__",     &PyFiltration::operator[])
        .def("index",           &PyFiltration::index)
        .def("__contains__",    &PyFiltration::contains)
        .def("__iter__",        [](const PyFiltration& f) { return py::make_iterator(f.begin(), f.end()); },
                                py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */)
        .def("sort",            [data_dim_cmp](PyFiltration& f) { f.sort(data_dim_cmp); })
        .def("sort",            [](PyFiltration& f, std::function<int (const PySimplex&, const PySimplex&)> cmp)
                                { f.sort([cmp](const PySimplex& s1, const PySimplex& s2) { return cmp(s1, s2) < 0; }); })
        .def("__repr__",        [](const PyFiltration& f) { std::ostringstream oss; oss << "Filtration with " << f.size() << " simplices"; return oss.str(); })
    ;
}
