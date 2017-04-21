#include <pybind11/pybind11.h>
namespace py = pybind11;

void init_simplex(py::module&);
void init_filtration(py::module&);
void init_rips(py::module&);
void init_freudenthal(py::module&);

void init_field(py::module&);
void init_persistence(py::module&);
void init_omnifield_persistence(py::module&);
void init_cohomology_persistence(py::module&);
void init_zigzag_persistence(py::module&);
void init_diagram(py::module&);
void init_distances(py::module&);

PYBIND11_PLUGIN(_dionysus)
{
    py::module m("_dionysus", "Dionysus python bindings");

    init_simplex(m);
    init_filtration(m);
    init_rips(m);
    init_freudenthal(m);

    init_field(m);
    init_persistence(m);
    init_cohomology_persistence(m);
    init_omnifield_persistence(m);
    init_zigzag_persistence(m);
    init_diagram(m);
    init_distances(m);

    return m.ptr();
}
