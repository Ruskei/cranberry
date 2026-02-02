#pragma once

#include "runtime_field.hpp"

struct Level {
  RuntimeField guess, target, residual, jacobi;
  const double scale_factor;
};
