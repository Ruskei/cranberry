#pragma once

#include "fdtd_types.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
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

struct RuntimeField {
  std::vector<double> v;
  int sx;
  int sy;
  int sz;
  RuntimeField(int sx, int sy, int sz);

  double &operator()(int x, int y, int z);
  const double &operator()(int x, int y, int z) const;
  RuntimeField &operator+=(const RuntimeField &other);
  // this = this + a * other
  void add_multiplied(double a, const RuntimeField &other);
  // this = other + a * this
  void multiply_add(double a, const RuntimeField &other);
  double dot(const RuntimeField &other) const;
  double norm2() const;

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
                            double scale_factor, double omega = 2.0 / 3.0);

void calculate_residual(const RuntimeField &guess, const RuntimeField &target,
                        RuntimeField &residual, double scale_factor);

void restrict(const RuntimeField &fine, RuntimeField &coarse);

void prolong(const RuntimeField &coarse, RuntimeField &fine);

bool reached_coarsest(int sx, int sy, int sz);

void solve_multigrid(MultigridContext &context, size_t level_idx);

MultigridContext create_multigrid_context(int sx, int sy, int sz,
                                          int smoothing);

void calculate_residual(const RuntimeField &guess, const RuntimeField &target,
                        RuntimeField &residual);

void update_residual_cg(RuntimeField &residual, double alpha,
                        const RuntimeField &search);

void apply_laplacian(const RuntimeField &guess, RuntimeField &out);

template <int N>
void solve_potential_cg(const Field<N, Component::Charge> &charge,
                        Field<N, Component::Potential> &potential) {
  const double tol = 1e-11;

  const int sx = charge.nx();
  const int sy = charge.ny();
  const int sz = charge.nz();

  const int n = sx * sy * sz;

  RuntimeField target{sx, sy, sz};
  target.read_from(charge);

  RuntimeField guess{sx, sy, sz};
  RuntimeField residual{sx, sy, sz};
  calculate_residual(guess, target, residual);

  RuntimeField search = residual; // copy
  RuntimeField laplacian_search{sx, sy, sz};

  double residual_norm = residual.norm2();

  for (auto k{0}; k < n; ++k) {
    apply_laplacian(search, laplacian_search);
    const double denom = search.dot(laplacian_search);
    if (std::abs(denom) < tol)
      break;

    const double alpha = residual_norm / denom;
    guess.add_multiplied(alpha, search);
    residual.add_multiplied(-alpha, laplacian_search);

    const double residual_norm_new = residual.norm2();
    if (std::abs(residual_norm_new) < tol)
      break;

    const double beta = residual_norm_new / residual_norm;
    search.multiply_add(beta, residual);

    residual_norm = residual_norm_new;
  }

  guess.write_into(potential);
}

template <int N>
void solve_potential_pcg(const Field<N, Component::Charge> &charge,
                         Field<N, Component::Potential> &potential) {
  const double tiny = 1e-30;
  const int smoothing = 3;
  const double tol = 1e-11;

  const int sx = charge.nx();
  const int sy = charge.ny();
  const int sz = charge.nz();

  assert(sx % 2 == 0);
  assert(sy % 2 == 0);
  assert(sz % 2 == 0);

  const int n = sx * sy * sz;

  RuntimeField target{sx, sy, sz};
  target.read_from(charge);

  RuntimeField guess{sx, sy, sz};
  RuntimeField residual{sx, sy, sz};

  calculate_residual(guess, target, residual);

  double target_norm2 = target.norm2();
  if (target_norm2 == 0.0)
    target_norm2 = 1.0;

  double residual_norm2 = residual.norm2();
  if (residual_norm2 <= (tol * tol) * target_norm2)
    return;

  RuntimeField conditioned{sx, sy, sz};

  MultigridContext context = create_multigrid_context(sx, sy, sz, smoothing);

  Level &first_level = context.levels[0];
  first_level.target = residual;
  solve_multigrid(context, 0);
  conditioned = first_level.guess;

  RuntimeField search{sx, sy, sz};
  RuntimeField laplacian_search{sx, sy, sz};
  search = conditioned;

  double rz_old = residual.dot(conditioned);

  int iterations = 0;
  for (auto k{0}; k < n; ++k) {
    iterations++;
    apply_laplacian(search, laplacian_search);
    const double denom = search.dot(laplacian_search);
    if (std::abs(denom) < tiny)
      break;

    const double alpha = rz_old / denom;

    guess.add_multiplied(alpha, search);
    residual.add_multiplied(-alpha, laplacian_search);

    if (residual.norm2() <= (tol * tol) * target_norm2)
      break;

    std::fill(first_level.guess.v.begin(), first_level.guess.v.end(), 0.0);
    first_level.target = residual;
    solve_multigrid(context, 0);
    conditioned = first_level.guess;

    const double rz_new = residual.dot(conditioned);

    if (std::abs(rz_old) < tiny)
      break;
    const double beta = rz_new / rz_old;
    search.multiply_add(beta, conditioned);

    rz_old = rz_new;
  }

  guess.write_into(potential);
  std::cout << "took " << iterations << " iterations!" << std::endl;
}
