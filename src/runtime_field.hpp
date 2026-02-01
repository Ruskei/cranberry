#pragma once

#include <assert.h>

#include "fdtd_types.hpp"

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
  // this = this + a * other
  void add_multiplied(double a, const RuntimeField &other) {
    assert(v.size() == other.v.size());
    for (size_t i{0}; i < v.size(); ++i)
      v[i] += a * other.v[i];
  }
  // this = other + a * this
  void multiply_add(double a, const RuntimeField &other) {
    assert(v.size() == other.v.size());
    for (size_t i{0}; i < v.size(); ++i)
      v[i] = other.v[i] + v[i] * a;
  }
  double dot(const RuntimeField &other) const {
    assert(v.size() == other.v.size());
    double sum = 0.0;
    for (size_t i{0}; i < v.size(); ++i)
      sum += v[i] * other.v[i];

    return sum;
  }
  double norm2() const {
    double sum = 0.0;
    for (size_t i{0}; i < v.size(); ++i)
      sum += v[i] * v[i];

    return sum;
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
