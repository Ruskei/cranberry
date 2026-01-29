#pragma once

#include <cmath>
#include <vector>

// J = E here but we keep duplicate for clarity
enum class Component { Ex, Ey, Ez, Jx, Jy, Jz, Hx, Hy, Hz, Charge, Potential };

template <int N, Component C> struct Extents;

template <int N> struct Extents<N, Component::Ex> {
  static constexpr int nx = N - 1, ny = N, nz = N;
  static constexpr double ox = 0.5, oy = 0.0, oz = 0.0;
};
template <int N> struct Extents<N, Component::Ey> {
  static constexpr int nx = N, ny = N - 1, nz = N;
  static constexpr double ox = 0.0, oy = 0.5, oz = 0.0;
};
template <int N> struct Extents<N, Component::Ez> {
  static constexpr int nx = N, ny = N, nz = N - 1;
  static constexpr double ox = 0.0, oy = 0.0, oz = 0.5;
};

template <int N> struct Extents<N, Component::Jx> {
  static constexpr int nx = N - 1, ny = N, nz = N;
  static constexpr double ox = 0.5, oy = 0.0, oz = 0.0;
};
template <int N> struct Extents<N, Component::Jy> {
  static constexpr int nx = N, ny = N - 1, nz = N;
  static constexpr double ox = 0.0, oy = 0.5, oz = 0.0;
};
template <int N> struct Extents<N, Component::Jz> {
  static constexpr int nx = N, ny = N, nz = N - 1;
  static constexpr double ox = 0.0, oy = 0.0, oz = 0.5;
};

template <int N> struct Extents<N, Component::Hx> {
  static constexpr int nx = N, ny = N - 1, nz = N - 1;
  static constexpr double ox = 0.0, oy = 0.5, oz = 0.5;
};
template <int N> struct Extents<N, Component::Hy> {
  static constexpr int nx = N - 1, ny = N, nz = N - 1;
  static constexpr double ox = 0.5, oy = 0.0, oz = 0.5;
};
template <int N> struct Extents<N, Component::Hz> {
  static constexpr int nx = N - 1, ny = N - 1, nz = N;
  static constexpr double ox = 0.5, oy = 0.5, oz = 0.0;
};

template <int N> struct Extents<N, Component::Charge> {
  static constexpr int nx = N - 1, ny = N - 1, nz = N - 1;
  static constexpr double ox = 0.0, oy = 0.0, oz = 0.0;
};
template <int N> struct Extents<N, Component::Potential> {
  static constexpr int nx = N - 1, ny = N - 1, nz = N - 1;
  static constexpr double ox = 0.0, oy = 0.0, oz = 0.0;
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

  double interpolate_at(double x, double y, double z) const {
    using e = Extents<N, C>;
    const double xp = x - e::ox;
    const double yp = y - e::oy;
    const double zp = z - e::oz;

    const int i = std::floor(xp);
    const int j = std::floor(yp);
    const int k = std::floor(zp);

    const double ax = xp - i, ay = yp - j, az = zp - k;

    const double wx[2] = {1.0 - ax, ax};
    const double wy[2] = {1.0 - ay, ay};
    const double wz[2] = {1.0 - az, az};

    const int sa = std::max(0, i) - i;
    const int sb = std::max(0, j) - j;
    const int sc = std::max(0, k) - k;

    const int ma = std::min(i + 2, nx()) - i;
    const int mb = std::min(j + 2, ny()) - j;
    const int mc = std::min(k + 2, nz()) - k;

    double sum = 0.0;
    double w_sum = 0.0;
    for (int a = sa; a < ma; ++a)
      for (int b = sb; b < mb; ++b)
        for (int c = sc; c < mc; ++c) {
          const double w = wx[a] * wy[b] * wz[c];
          w_sum += w;
          sum += w * (*this)(i + a, j + b, k + c);
        }

    if (w_sum > 1e-11)
      sum /= w_sum;

    return sum;
  }
};

template <int N> struct EField {
  Field<N, Component::Ex> x;
  Field<N, Component::Ex> cxe;
  Field<N, Component::Ex> cxh;

  Field<N, Component::Ey> y;
  Field<N, Component::Ey> cye;
  Field<N, Component::Ey> cyh;

  Field<N, Component::Ez> z;
  Field<N, Component::Ez> cze;
  Field<N, Component::Ez> czh;
};

template <int N> struct HField {
  Field<N, Component::Hx> x;
  Field<N, Component::Hx> cxh;
  Field<N, Component::Hx> cxe;

  Field<N, Component::Hy> y;
  Field<N, Component::Hy> cyh;
  Field<N, Component::Hy> cye;

  Field<N, Component::Hz> z;
  Field<N, Component::Hz> czh;
  Field<N, Component::Hz> cze;
};

template <int N> struct JField {
  Field<N, Component::Jx> x;
  Field<N, Component::Jy> y;
  Field<N, Component::Jz> z;
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
