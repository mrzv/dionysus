#pragma once

#include <dionysus/diagram.h>

#include "simplex.h"

using PyIndex   = unsigned;
using PyDiagram = dionysus::Diagram<PySimplex::Data, PyIndex>;
