#include <sstream>

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "field.h"

void init_field(py::module& m)
{
    using namespace pybind11::literals;
    py::class_<PyZpField>(m, "Zp", "arithmetic mod p")
        .def(py::init<PyZpField::Element>())
        .def("__repr__",        [](const PyZpField& f) { std::ostringstream oss; oss << "Z mod " << f.prime(); return oss.str(); })
        .def("id",              &PyZpField::id,                      "1")
        .def("zero",            &PyZpField::zero,                    "0")
        .def("init",            &PyZpField::init,     "a"_a,         "(a % p + p) % p")
        .def("neg",             &PyZpField::neg,      "a"_a,         "-a")
        .def("add",             &PyZpField::add,      "a"_a, "b"_a,  "a + b")
        .def("inv",             &PyZpField::inv,      "a"_a,         "1/a")
        .def("mul",             &PyZpField::mul,      "a"_a, "b"_a,  "a * b")
        .def("div",             &PyZpField::div,      "a"_a, "b"_a,  "a / b")
        .def("is_zero",         &PyZpField::is_zero,  "a"_a,         "a == 0")
        .def("prime",           &PyZpField::prime,                    "p")
    ;

    py::class_<PyQElement>(m, "QElement", "rational number")
        .def_readonly("numerator",    &PyQElement::numerator,    "numerator")
        .def_readonly("denominator",  &PyQElement::denominator,  "denominator")
        .def("__repr__",        [](const PyQElement& e) { std::ostringstream oss; oss << e; return oss.str(); })
    ;
}
