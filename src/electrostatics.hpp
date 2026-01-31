#pragma once

#include "fdtd_types.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

template <int N> struct ResidualStats {
  double l2_abs;
  double l2_rel;
  double linf_abs;
};

template <int N>
ResidualStats<N>
calculate_residuals(const Field<N, Component::Charge> &charge,
                    const Field<N, Component::Potential> &potential) {
  double sum_r2 = 0.0;
  double sum_b2 = 0.0;
  double linf = 0.0;
  std::size_t count = 0;

  for (int x = 1; x < N - 2; ++x)
    for (int y = 1; y < N - 2; ++y)
      for (int z = 1; z < N - 2; ++z) {
        const double lap = potential(x - 1, y, z) + potential(x + 1, y, z) +
                           potential(x, y - 1, z) + potential(x, y + 1, z) +
                           potential(x, y, z - 1) + potential(x, y, z + 1) -
                           6.0 * potential(x, y, z);

        const double r = lap + charge(x, y, z);
        const double b = -charge(x, y, z);

        sum_r2 += r * r;
        sum_b2 += b * b;
        linf = std::max(linf, std::abs(r));
        ++count;
      }

  const double inv_n = (count > 0) ? (1.0 / static_cast<double>(count)) : 0.0;
  const double l2_abs = std::sqrt(sum_r2 * inv_n);

  const double denom = std::sqrt(sum_b2 * inv_n);
  const double l2_rel = (denom > std::numeric_limits<double>::epsilon())
                            ? (l2_abs / denom)
                            : l2_abs;

  return {l2_abs, l2_rel, linf};
}

template <int N>
void solve_potential(const Field<N, Component::Charge> &charge,
                     Field<N, Component::Potential> &potential) {
  for (auto i{0}; i < 500; ++i)
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

  RuntimeField &operator+=(const RuntimeField &other) {
    assert(sx == other.sx && sy == other.sy && sz == other.sz);
    for (size_t i{0}; i < v.size(); ++i)
      v[i] += other.v[i];

    return *this;
  }

  template <class Field> void write_into(Field &field) const {
    assert(sx == field.nx());
    assert(sy == field.ny());
    assert(sz == field.nz());
    for (size_t i{0}; i < v.size(); ++i)
      field.v[i] = v[i];
  }

  template <class Field> void read_from(const Field &field) {
    assert(sx == field.nx());
    assert(sy == field.ny());
    assert(sz == field.nz());
    for (size_t i{0}; i < v.size(); ++i)
      v[i] = field.v[i];
  }
};

struct Level {
  RuntimeField guess, target, residual;
  const double scale_factor;
};

struct MultigridContext {
  std::vector<Level> levels;
  const int smoothing{1};
};

void smooth_weighted_jacobi(RuntimeField &guess, const RuntimeField &target,
                            double scale_factor, double omega = 2.0 / 3.0) {
  assert(guess.sx == target.sx && guess.sy == target.sy &&
         guess.sz == target.sz);
  const int sx = guess.sx, sy = guess.sy, sz = guess.sz;
  const int yz = sy * sz;
  const double h2 = scale_factor * scale_factor;

  RuntimeField new_guess{sx, sy, sz};

  for (auto x{1}; x < sx - 1; ++x) {
    const int xb = x * yz;
    for (auto y{1}; y < sy - 1; ++y) {
      const int yb = xb + y * sz;
      for (auto z{1}; z < sz - 1; ++z) {
        const int i = yb + z;

        const double jacobi =
            (guess.v[i - yz] + guess.v[i + yz] + guess.v[i - sz] +
             guess.v[i + sz] + guess.v[i - 1] + guess.v[i + 1] +
             h2 * target.v[i]) /
            6.0;

        new_guess.v[i] = (1.0 - omega) * guess.v[i] + omega * jacobi;
      }
    }
  }

  guess.v.swap(new_guess.v);
}

void calculate_residual(const RuntimeField &guess, const RuntimeField &target,
                        RuntimeField &residual, double scale_factor) {
  const int sx = guess.sx, sy = guess.sy, sz = guess.sz;
  const int yz = sy * sz;
  const double c = 1.0 / (scale_factor * scale_factor);

  for (auto x{1}; x < sx - 1; ++x) {
    const int xb = x * yz;
    for (auto y{1}; y < sy - 1; ++y) {
      const int yb = xb + y * sz;
      for (auto z{1}; z < sz - 1; ++z) {
        const int i = yb + z;
        residual.v[i] = target.v[i] - c * (6.0 * guess.v[i] -
                                           (guess.v[i - yz] + guess.v[i + yz] +
                                            guess.v[i - sz] + guess.v[i + sz] +
                                            guess.v[i - 1] + guess.v[i + 1]));
      }
    }
  }
}

void restrict(const RuntimeField &fine, RuntimeField &coarse) {
  assert(fine.sx % 2 == 0);
  assert(fine.sy % 2 == 0);
  assert(fine.sz % 2 == 0);

  auto w = [](int d) -> double { return (d == 0) ? 2.0 : 1.0; };

  for (int i = 1; i < coarse.sx - 1; ++i)
    for (int j = 1; j < coarse.sy - 1; ++j)
      for (int k = 1; k < coarse.sz - 1; ++k) {
        const int i0 = 2 * i, j0 = 2 * j, k0 = 2 * k;

        assert(i0 + 1 < fine.sx && j0 + 1 < fine.sy && k0 + 1 < fine.sz);

        double sum = 0.0;
        for (int a = -1; a <= 1; ++a)
          for (int b = -1; b <= 1; ++b)
            for (int c = -1; c <= 1; ++c)
              sum += w(a) * w(b) * w(c) * fine(i0 + a, j0 + b, k0 + c);

        coarse(i, j, k) = sum / 64.0;
      }
}

void prolong(const RuntimeField &coarse, RuntimeField &fine) {
  for (int i = 1; i < fine.sx - 1; ++i) {
    const int I = i / 2;
    const int di = i % 2;
    const int I1 = std::min(I + di, coarse.sx - 1);
    const double wx0 = di ? 0.5 : 1.0;
    const double wx1 = 1.0 - wx0;

    for (int j = 1; j < fine.sy - 1; ++j) {
      const int J = j / 2;
      const int dj = j % 2;
      const int J1 = std::min(J + dj, coarse.sy - 1);
      const double wy0 = dj ? 0.5 : 1.0;
      const double wy1 = 1.0 - wy0;

      for (int k = 1; k < fine.sz - 1; ++k) {
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
}

bool reached_coarsest(int sx, int sy, int sz) {
  return sx < 4 || sy < 4 || sz < 4;
}

void solve_multigrid(MultigridContext &context, size_t level_idx) {
  assert(level_idx < context.levels.size());
  Level &level = context.levels[level_idx];
  if (level_idx == context.levels.size() - 1) {
    for (auto i{0}; i < 50; ++i)
      smooth_weighted_jacobi(level.guess, level.target, level.scale_factor);

    return;
  }

  for (auto i{0}; i < context.smoothing; ++i)
    smooth_weighted_jacobi(level.guess, level.target, level.scale_factor);

  Level &next_level = context.levels[level_idx + 1];
  std::fill(next_level.guess.v.begin(), next_level.guess.v.end(), 0);
  calculate_residual(level.guess, level.target, level.residual,
                     level.scale_factor);
  restrict(level.residual, next_level.target);
  solve_multigrid(context, level_idx + 1);

  prolong(next_level.guess, level.residual);
  level.guess += level.residual;

  for (auto i{0}; i < context.smoothing; ++i)
    smooth_weighted_jacobi(level.guess, level.target, level.scale_factor);
}

MultigridContext create_multigrid_context(int sx, int sy, int sz,
                                          int smoothing) {
  std::vector<Level> levels;
  double scale_factor = 1.0;
  while (!reached_coarsest(sx, sy, sz)) {
    levels.push_back(Level{
        RuntimeField{sx, sy, sz},
        RuntimeField{sx, sy, sz},
        RuntimeField{sx, sy, sz},
        scale_factor,
    });
    sx /= 2;
    sy /= 2;
    sz /= 2;
    scale_factor *= 2.0;
  }

  return MultigridContext{std::move(levels), smoothing};
}

template <int N>
void solve_potential_multigrid(const Field<N, Component::Charge> &charge,
                               Field<N, Component::Potential> &potential,
                               int iterations) {
  const int smoothing = 3;

  const int sx = charge.nx();
  const int sy = charge.ny();
  const int sz = charge.nz();

  assert(sx % 2 == 0);
  assert(sy % 2 == 0);
  assert(sz % 2 == 0);

  MultigridContext context = create_multigrid_context(sx, sy, sz, smoothing);
  Level &first_level = context.levels[0];
  first_level.target.read_from(charge);
  for (auto i{0}; i < iterations; ++i)
    solve_multigrid(context, 0);
  first_level.guess.write_into(potential);
}
