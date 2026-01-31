#pragma once

#include "fdtd_types.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

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
  for (auto i{0}; i < 50; ++i)
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

struct RuntimeField {
  std::vector<double> v;
  int sx;
  int sy;
  int sz;
  RuntimeField(int sx, int sy, int sz)
      : v{std::vector<double>(sx * sy * sz)}, sx{sx}, sy{sy}, sz{sz} {}

  double &operator()(int x, int y, int z) {
    return v[x * sy * sz + y * sz + z];
  }
  const double &operator()(int x, int y, int z) const {
    return v[x * sy * sz + y * sz + z];
  }

  RuntimeField operator-(const RuntimeField &other) const {
    assert(sx == other.sx && sy == other.sy && sz == other.sz);
    RuntimeField result = RuntimeField{sx, sy, sz};
    for (size_t i{0}; i < v.size(); ++i)
      result.v[i] = v[i] - other.v[i];

    return result;
  }

  RuntimeField &operator+=(const RuntimeField &other) {
    assert(sx == other.sx && sy == other.sy && sz == other.sz);
    for (size_t i{0}; i < v.size(); ++i)
      v[i] += other.v[i];

    return *this;
  }

  template <class Field> void write_into(Field &field) {
    assert(sx == field.nx());
    assert(sy == field.ny());
    assert(sz == field.nz());
    for (size_t i{0}; i < v.size(); ++i)
      field.v[i] = v[i];
  }
};

template <class Field> RuntimeField runtime_field_from(Field &field) {
  RuntimeField runtime_field = RuntimeField{field.nx(), field.ny(), field.nz()};
  for (size_t i{0}; i < runtime_field.v.size(); ++i)
    runtime_field.v[i] = field.v[i];

  return runtime_field;
}

// smooths -∇^2 guess = target
void smooth(RuntimeField &guess, const RuntimeField &target,
            double scale_factor) {
  assert(guess.sx == target.sx && guess.sy == target.sy &&
         guess.sz == target.sz);
  const double c = scale_factor * scale_factor;
  for (auto x{1}; x < guess.sx - 1; ++x)
    for (auto y{1}; y < guess.sy - 1; ++y)
      for (auto z{1}; z < guess.sz - 1; ++z) {
        guess(x, y, z) =
            (guess(x - 1, y, z) + guess(x + 1, y, z) + guess(x, y - 1, z) +
             guess(x, y + 1, z) + guess(x, y, z - 1) + guess(x, y, z + 1) +
             c * target(x, y, z)) /
            6.0;
      }
}

RuntimeField laplacian(RuntimeField &field, double scale_factor) {
  RuntimeField result = RuntimeField(field.sx, field.sy, field.sz);
  const double c = 1.0 / (scale_factor * scale_factor);

  for (auto x{1}; x < result.sx - 1; ++x)
    for (auto y{1}; y < result.sy - 1; ++y)
      for (auto z{1}; z < result.sz - 1; ++z) {
        result(x, y, z) = c * (6.0 * field(x, y, z) -
                               (field(x - 1, y, z) + field(x + 1, y, z) +
                                field(x, y - 1, z) + field(x, y + 1, z) +
                                field(x, y, z - 1) + field(x, y, z + 1)));
      }

  return result;
}

RuntimeField restrict(const RuntimeField &fine) {
  assert(fine.sx % 2 == 0);
  assert(fine.sy % 2 == 0);
  assert(fine.sz % 2 == 0);

  RuntimeField coarse{fine.sx / 2, fine.sy / 2, fine.sz / 2};

  auto w = [](int d) -> double { return (d == 0) ? 2.0 : 1.0; };

  for (int i = 0; i < coarse.sx; ++i)
    for (int j = 0; j < coarse.sy; ++j)
      for (int k = 0; k < coarse.sz; ++k) {
        const int i0 = 2 * i, j0 = 2 * j, k0 = 2 * k;

        if (i0 - 1 < 0 || i0 + 1 >= fine.sx || j0 - 1 < 0 ||
            j0 + 1 >= fine.sy || k0 - 1 < 0 || k0 + 1 >= fine.sz) {
          coarse(i, j, k) = fine(i0, j0, k0);
          continue;
        }

        double sum = 0.0;
        for (int a = -1; a <= 1; ++a)
          for (int b = -1; b <= 1; ++b)
            for (int c = -1; c <= 1; ++c)
              sum += w(a) * w(b) * w(c) * fine(i0 + a, j0 + b, k0 + c);

        coarse(i, j, k) = sum / 64.0;
      }

  return coarse;
}

RuntimeField prolong(const RuntimeField &coarse) {
  RuntimeField fine{coarse.sx * 2, coarse.sy * 2, coarse.sz * 2};

  for (int i = 0; i < fine.sx; ++i) {
    const int I = i / 2;
    const int di = i % 2;
    const int I1 = std::min(I + di, coarse.sx - 1);
    const double wx0 = di ? 0.5 : 1.0;
    const double wx1 = 1.0 - wx0;

    for (int j = 0; j < fine.sy; ++j) {
      const int J = j / 2;
      const int dj = j % 2;
      const int J1 = std::min(J + dj, coarse.sy - 1);
      const double wy0 = dj ? 0.5 : 1.0;
      const double wy1 = 1.0 - wy0;

      for (int k = 0; k < fine.sz; ++k) {
        const int K = k / 2;
        const int dk = k % 2;
        const int K1 = std::min(K + dk, coarse.sz - 1);
        const double wz0 = dk ? 0.5 : 1.0;
        const double wz1 = 1.0 - wz0;

        fine(i, j, k) = wx0 * wy0 * wz0 * coarse(I, J, K) +
                        wx1 * wy0 * wz0 * coarse(I1, J, K) +
                        wx0 * wy1 * wz0 * coarse(I, J1, K) +
                        wx1 * wy1 * wz0 * coarse(I1, J1, K) +
                        wx0 * wy0 * wz1 * coarse(I, J, K1) +
                        wx1 * wy0 * wz1 * coarse(I1, J, K1) +
                        wx0 * wy1 * wz1 * coarse(I, J1, K1) +
                        wx1 * wy1 * wz1 * coarse(I1, J1, K1);
      }
    }
  }

  return fine;
}

// a scale_factor=2 means the grid has half the cells, meaning each cell is
// twice as big
RuntimeField solve_multigrid(RuntimeField guess, const RuntimeField &target,
                             double scale_factor, int smoothing) {
  assert(target.sx == guess.sx && target.sy == guess.sy &&
         target.sz == guess.sz);

  if (guess.sx < 4 || guess.sy < 4 || guess.sz < 4) {
    for (auto i{0}; i < 5; ++i)
      smooth(guess, target, scale_factor);

    return guess;
  }

  for (auto i{0}; i < smoothing; ++i)
    smooth(guess, target, scale_factor);

  RuntimeField residual = restrict(target - laplacian(guess, scale_factor));

  RuntimeField correction =
      solve_multigrid(RuntimeField{residual.sx, residual.sy, residual.sz},
                      residual, scale_factor * 2.0, smoothing);
  // w-cycle
  // correction = solve_multigrid(correction, residual, scale_factor * 2.0, smoothing);

  guess += prolong(correction);

  for (auto i{0}; i < smoothing; ++i)
    smooth(guess, target, scale_factor);

  return guess;
}

template <int N>
void solve_potential_multigrid(const Field<N, Component::Charge> &charge,
                               Field<N, Component::Potential> &potential) {
  const int smoothing = 10;
  RuntimeField charge_runtime = runtime_field_from(charge);
  const int sx = charge_runtime.sx;
  const int sy = charge_runtime.sy;
  const int sz = charge_runtime.sz;

  assert(sx % 2 == 0);
  assert(sy % 2 == 0);
  assert(sz % 2 == 0);

  RuntimeField potential_runtime =
      solve_multigrid(RuntimeField{sx, sy, sz}, charge_runtime, 1.0, smoothing);
  potential_runtime.write_into(potential);
}
