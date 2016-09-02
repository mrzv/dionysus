#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "simplex.h"
#include "filtration.h"

void init_filtration(py::module& m)
{
    py::class_<PyFiltration>(m, "Filtration")
        .def(py::init<>())
        .def("add", [](PyFiltration* f, const PySimplex& s) { f->push_back(s); });
    ;
}
