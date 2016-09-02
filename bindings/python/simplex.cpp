#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "simplex.h"

void init_simplex(py::module& m)
{
    py::class_<PySimplex>(m, "Simplex")
        .def(py::init<>())
    ;
}
