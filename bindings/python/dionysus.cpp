#include <pybind11/pybind11.h>
namespace py = pybind11;

void init_simplex(py::module&);
void init_filtration(py::module&);

PYBIND11_PLUGIN(_dionysus)
{
    py::module m("_dionysus", "Dionysus python bindings");

    init_simplex(m);
    init_filtration(m);

    return m.ptr();
}
