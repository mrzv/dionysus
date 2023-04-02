#pragma once

#include <dionysus/fields/zp.h>
#include <dionysus/fields/q.h>

using PyZpField = dionysus::ZpField<long>;      // long to be usable with Q in OmniFieldPersistence
using PyQElement = dionysus::Q<long>::Element;
