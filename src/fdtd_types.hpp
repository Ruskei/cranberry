#pragma once

#include <vector>

enum class Component { Ex, Ey, Ez, Hx, Hy, Hz };

template <int N, Component C> struct Extents;

template <int N> struct Extents<N, Component::Ex> {
  static constexpr int nx = N - 1, ny = N, nz = N;
};
template <int N> struct Extents<N, Component::Ey> {
  static constexpr int nx = N, ny = N - 1, nz = N;
};
template <int N> struct Extents<N, Component::Ez> {
  static constexpr int nx = N, ny = N, nz = N - 1;
};

template <int N> struct Extents<N, Component::Hx> {
  static constexpr int nx = N, ny = N - 1, nz = N - 1;
};
template <int N> struct Extents<N, Component::Hy> {
  static constexpr int nx = N - 1, ny = N, nz = N - 1;
};
template <int N> struct Extents<N, Component::Hz> {
  static constexpr int nx = N - 1, ny = N - 1, nz = N;
};

template <int N, Component C> struct Layout {
  static constexpr int nx = Extents<N, C>::nx;
  static constexpr int ny = Extents<N, C>::ny;
  static constexpr int nz = Extents<N, C>::nz;

  static constexpr int sx = ny * nz;
  static constexpr int sy = nz;
  static constexpr int sz = 1;

  static constexpr int size = nx * ny * nz;

  static constexpr int idx(int x, int y, int z) { return x * sx + y * sy + z; }
};

template <int N, Component C> struct Field {
  std::vector<double> v;
  Field() : v(Layout<N, C>::size) {}
  double &operator()(int x, int y, int z) {
    return v[Layout<N, C>::idx(x, y, z)];
  }
  const double &operator()(int x, int y, int z) const {
    return v[Layout<N, C>::idx(x, y, z)];
  }

  static constexpr int nx() { return Layout<N, C>::nx; }
  static constexpr int ny() { return Layout<N, C>::ny; }
  static constexpr int nz() { return Layout<N, C>::nz; }
};

template <int N, Component C>
void setup_field_coeffients(Field<N, C> &ca, Field<N, C> &cb, double factor) {
  using l = Layout<N, C>;
  for (auto x{0}; x < l::nx; ++x)
    for (auto y{0}; y < l::ny; ++y)
      for (auto z{0}; z < l::nz; ++z) {
        ca(x, y, z) = 1.0;
        cb(x, y, z) = factor;
      }
}
