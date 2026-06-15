#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <dionysus/boundary-matrix.h>

#include "filtration.h"
#include "persistence.h"

template<class PyFiltration>
PyMatrixFiltration boundary(const PyFiltration& f)
{
    return dionysus::make_boundary_matrix_filtration<PyMatrixFiltration>(f);
}

template<class PyFiltration>
PyMatrixFiltration coboundary(const PyFiltration& f)
{
    return dionysus::make_coboundary_matrix_filtration<PyMatrixFiltration>(f);
}


void init_boundary(py::module& m)
{
    m.def("boundary", &boundary<PyFiltration>, "compute boundary matrix of the filtration");
    m.def("coboundary", &coboundary<PyFiltration>, "compute coboundary matrix of the filtration");
    m.def("boundary", &boundary<PyMultiFiltration>, "compute boundary matrix of the filtration");
    m.def("coboundary", &coboundary<PyMultiFiltration>, "compute coboundary matrix of the filtration");
}
