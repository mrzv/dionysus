#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <stdexcept>
namespace py = pybind11;

#include <dionysus/freudenthal.h>

#include "filtration.h"

template<class T>
PyFiltration fill_freudenthal_(py::array a, bool reverse)
{
    std::vector<int> shape(a.shape(), a.shape() + a.ndim());
    std::vector<size_t> strides;
    strides.reserve(a.ndim());
    for (size_t i = 0; i < a.ndim(); ++i)
    {
        if (a.strides(i) < 0)
            throw std::runtime_error("negative array strides are not supported");
        strides.push_back(a.strides(i) / sizeof(T));
    }

    const T* a_data = static_cast<const T*>(a.data(0));
    return dionysus::fill_freudenthal<PyFiltration>(shape, strides, a_data, reverse, DataDimCmp(reverse));
}

PyFiltration fill_freudenthal(py::array a, bool reverse)
{
    if (a.dtype().is(py::dtype::of<float>()))
        return fill_freudenthal_<float>(a, reverse);
    else if (a.dtype().is(py::dtype::of<double>()))
        return fill_freudenthal_<double>(a, reverse);
    else
        throw std::runtime_error("Unknown array dtype");
}

void init_freudenthal(py::module& m)
{
    using namespace pybind11::literals;
    m.def("fill_freudenthal",  &fill_freudenthal,
          "data"_a, "reverse"_a = false,
          "returns (sorted) lower-star (or upper-star if ``reverse = True``) filtration filled with the Freudenthal triangulation of the grid in the array `data`");
}
