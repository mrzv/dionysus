#include <vector>

#include <dionysus/reduced-matrix.h>
#include <dionysus/matrix-filtration.h>
#include <dionysus/trails-chains.h>

#include "filtration.h"
#include "field.h"

using PyReducedMatrix = dionysus::ReducedMatrix<PyZpField>;
using PyReducedMatrixNoNegative = dionysus::OrdinaryPersistenceNoNegative<PyZpField>;
using PyReducedMatrixWithV = dionysus::OrdinaryPersistenceWithV<PyZpField>;
using PyReducedMatrixNoNegativeWithV = dionysus::OrdinaryPersistenceNoNegativeWithV<PyZpField>;

using PyMatrixFiltration = dionysus::MatrixFiltration<PyReducedMatrix,PySimplex::Data>;
using Dimensions = PyMatrixFiltration::Dimensions;
using Values = PyMatrixFiltration::Values;
