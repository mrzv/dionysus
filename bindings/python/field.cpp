#include <sstream>

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "field.h"

void init_field(py::module& m)
{
    py::class_<PyZpField>(m, "Zp", "arithmetic mod p")
        .def(py::init<PyZpField::Element>())
        .def("__repr__",        [](const PyZpField& f) { std::ostringstream oss; oss << "Z mod " << f.prime(); return oss.str(); })
    ;

    py::class_<PyQElement>(m, "QElement", "rational number")
        .def_readonly("numerator",    &PyQElement::numerator,    "numerator")
        .def_readonly("denominator",  &PyQElement::denominator,  "denominator")
        .def("__repr__",        [](const PyQElement& e) { std::ostringstream oss; oss << e; return oss.str(); })
    ;
}
