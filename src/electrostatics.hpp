#pragma once

#include "fdtd_types.hpp"

template <int N>
double calculate_residuals(const Field<N, Component::Charge> &charge,
                           const Field<N, Component::Potential> &potential) {
  double max_residual = 0.0;
  for (auto x{1}; x < N - 2; ++x)
    for (auto y{1}; y < N - 2; ++y)
      for (auto z{1}; z < N - 2; ++z) {
        const double lap = potential(x - 1, y, z) + potential(x + 1, y, z) +
                           potential(x, y - 1, z) + potential(x, y + 1, z) +
                           potential(x, y, z - 1) + potential(x, y, z + 1) -
                           6.0 * potential(x, y, z);
        const double residual = std::abs(lap + charge(x, y, z));
        max_residual = std::max(max_residual, residual);
      }

  return max_residual;
}

template <int N>
void solve_potential(const Field<N, Component::Charge> &charge,
                     Field<N, Component::Potential> &potential) {
  for (auto i{0}; i < 1'000; ++i)
    for (auto x{1}; x < N - 2; ++x)
      for (auto y{1}; y < N - 2; ++y)
        for (auto z{1}; z < N - 2; ++z) {
          potential(x, y, z) =
              (potential(x - 1, y, z) + potential(x + 1, y, z) +
               potential(x, y - 1, z) + potential(x, y + 1, z) +
               potential(x, y, z - 1) + potential(x, y, z + 1) +
               charge(x, y, z)) /
              6.0;
        }
}
