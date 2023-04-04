#include <vector>

#include <dionysus/reduced-matrix.h>
#include <dionysus/matrix-filtration.h>

#include "filtration.h"
#include "field.h"

using PyReducedMatrix = dionysus::ReducedMatrix<PyZpField>;
using PyMatrixFiltration = dionysus::MatrixFiltration<PyReducedMatrix>;
using Dimensions = std::vector<short unsigned>;
