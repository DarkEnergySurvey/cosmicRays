#pragma once
#define WRAP(x) void wrap ## x(py::module_ &mod)

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include <pybind11/stl_bind.h>
