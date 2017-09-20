#pragma once

#include <dionysus/filtration.h>

#include "simplex.h"

using PyFiltration = dionysus::Filtration<PySimplex, bmi::hashed_unique<bmi::identity<PySimplex>>, true>;
