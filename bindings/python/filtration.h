#pragma once

#include <dionysus/filtration.h>
#include <dionysus/multi-filtration.h>
#include <dionysus/linked-multi-filtration.h>

#include "simplex.h"

using PyFiltration = dionysus::Filtration<PySimplex, bmi::hashed_unique<bmi::identity<PySimplex>>, true>;
using PyMultiFiltration = dionysus::MultiFiltration<PySimplex, true>;
using PyLinkedMultiFiltration = dionysus::LinkedMultiFiltration<PySimplex, true>;
