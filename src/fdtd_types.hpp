#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

// J = E here but we keep duplicate for clarity
enum class Component { Ex, Ey, Ez, Jx, Jy, Jz, Hx, Hy, Hz, Charge, Potential };

namespace grid_detail {
constexpr bool is_power_of_two(int value) {
  return value > 0 && (value & (value - 1)) == 0;
}

constexpr bool is_valid_grid_extent(int value) {
  return value >= 5 && is_power_of_two(value - 1);
}
} // namespace grid_detail

template <int NX, int NY, int NZ> struct GridShape {
  static_assert(grid_detail::is_valid_grid_extent(NX),
                "NX must be 2^n + 1 and >= 5");
  static_assert(grid_detail::is_valid_grid_extent(NY),
                "NY must be 2^n + 1 and >= 5");
  static_assert(grid_detail::is_valid_grid_extent(NZ),
                "NZ must be 2^n + 1 and >= 5");

  static constexpr int nx = NX;
  static constexpr int ny = NY;
  static constexpr int nz = NZ;
};

template <int NX, int NY, int NZ, Component C> struct Extents;

template <int NX, int NY, int NZ> struct Extents<NX, NY, NZ, Component::Ex> {
  using Shape = GridShape<NX, NY, NZ>;
  static constexpr int nx = NX - 1, ny = NY, nz = NZ;
  static constexpr double ox = 0.5, oy = 0.0, oz = 0.0;
};
template <int NX, int NY, int NZ> struct Extents<NX, NY, NZ, Component::Ey> {
  using Shape = GridShape<NX, NY, NZ>;
  static constexpr int nx = NX, ny = NY - 1, nz = NZ;
  static constexpr double ox = 0.0, oy = 0.5, oz = 0.0;
};
template <int NX, int NY, int NZ> struct Extents<NX, NY, NZ, Component::Ez> {
  using Shape = GridShape<NX, NY, NZ>;
  static constexpr int nx = NX, ny = NY, nz = NZ - 1;
  static constexpr double ox = 0.0, oy = 0.0, oz = 0.5;
};

template <int NX, int NY, int NZ> struct Extents<NX, NY, NZ, Component::Jx> {
  using Shape = GridShape<NX, NY, NZ>;
  static constexpr int nx = NX - 1, ny = NY, nz = NZ;
  static constexpr double ox = 0.5, oy = 0.0, oz = 0.0;
};
template <int NX, int NY, int NZ> struct Extents<NX, NY, NZ, Component::Jy> {
  using Shape = GridShape<NX, NY, NZ>;
  static constexpr int nx = NX, ny = NY - 1, nz = NZ;
  static constexpr double ox = 0.0, oy = 0.5, oz = 0.0;
};
template <int NX, int NY, int NZ> struct Extents<NX, NY, NZ, Component::Jz> {
  using Shape = GridShape<NX, NY, NZ>;
  static constexpr int nx = NX, ny = NY, nz = NZ - 1;
  static constexpr double ox = 0.0, oy = 0.0, oz = 0.5;
};

template <int NX, int NY, int NZ> struct Extents<NX, NY, NZ, Component::Hx> {
  using Shape = GridShape<NX, NY, NZ>;
  static constexpr int nx = NX, ny = NY - 1, nz = NZ - 1;
  static constexpr double ox = 0.0, oy = 0.5, oz = 0.5;
};
template <int NX, int NY, int NZ> struct Extents<NX, NY, NZ, Component::Hy> {
  using Shape = GridShape<NX, NY, NZ>;
  static constexpr int nx = NX - 1, ny = NY, nz = NZ - 1;
  static constexpr double ox = 0.5, oy = 0.0, oz = 0.5;
};
template <int NX, int NY, int NZ> struct Extents<NX, NY, NZ, Component::Hz> {
  using Shape = GridShape<NX, NY, NZ>;
  static constexpr int nx = NX - 1, ny = NY - 1, nz = NZ;
  static constexpr double ox = 0.5, oy = 0.5, oz = 0.0;
};

template <int NX, int NY, int NZ>
struct Extents<NX, NY, NZ, Component::Charge> {
  using Shape = GridShape<NX, NY, NZ>;
  static constexpr int nx = NX - 1, ny = NY - 1, nz = NZ - 1;
  static constexpr double ox = 0.0, oy = 0.0, oz = 0.0;
};
template <int NX, int NY, int NZ>
struct Extents<NX, NY, NZ, Component::Potential> {
  using Shape = GridShape<NX, NY, NZ>;
  static constexpr int nx = NX - 1, ny = NY - 1, nz = NZ - 1;
  static constexpr double ox = 0.0, oy = 0.0, oz = 0.0;
};

template <int NX, int NY, int NZ, Component C> struct Layout {
  using Shape = GridShape<NX, NY, NZ>;
  static constexpr int nx = Extents<NX, NY, NZ, C>::nx;
  static constexpr int ny = Extents<NX, NY, NZ, C>::ny;
  static constexpr int nz = Extents<NX, NY, NZ, C>::nz;

  static constexpr int sx = ny * nz;
  static constexpr int sy = nz;
  static constexpr int sz = 1;

  static constexpr int size = nx * ny * nz;

  static constexpr int idx(int x, int y, int z) { return x * sx + y * sy + z; }
};

template <int NX, int NY, int NZ, Component C> struct Field {
  using Shape = GridShape<NX, NY, NZ>;
  std::vector<double> v;
  Field() : v(Layout<NX, NY, NZ, C>::size) {}
  double &operator()(int x, int y, int z) {
    return v[Layout<NX, NY, NZ, C>::idx(x, y, z)];
  }
  const double &operator()(int x, int y, int z) const {
    return v[Layout<NX, NY, NZ, C>::idx(x, y, z)];
  }

  static constexpr int nx() { return Layout<NX, NY, NZ, C>::nx; }
  static constexpr int ny() { return Layout<NX, NY, NZ, C>::ny; }
  static constexpr int nz() { return Layout<NX, NY, NZ, C>::nz; }

  double interpolate_at(double x, double y, double z) const {
    using e = Extents<NX, NY, NZ, C>;
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

template <int NX, int NY, int NZ> struct EField {
  using Shape = GridShape<NX, NY, NZ>;
  Field<NX, NY, NZ, Component::Ex> x;
  Field<NX, NY, NZ, Component::Ex> cxe;
  Field<NX, NY, NZ, Component::Ex> cxh;

  Field<NX, NY, NZ, Component::Ey> y;
  Field<NX, NY, NZ, Component::Ey> cye;
  Field<NX, NY, NZ, Component::Ey> cyh;

  Field<NX, NY, NZ, Component::Ez> z;
  Field<NX, NY, NZ, Component::Ez> cze;
  Field<NX, NY, NZ, Component::Ez> czh;
};

template <int NX, int NY, int NZ> struct HField {
  using Shape = GridShape<NX, NY, NZ>;
  Field<NX, NY, NZ, Component::Hx> x;
  Field<NX, NY, NZ, Component::Hx> cxh;
  Field<NX, NY, NZ, Component::Hx> cxe;

  Field<NX, NY, NZ, Component::Hy> y;
  Field<NX, NY, NZ, Component::Hy> cyh;
  Field<NX, NY, NZ, Component::Hy> cye;

  Field<NX, NY, NZ, Component::Hz> z;
  Field<NX, NY, NZ, Component::Hz> czh;
  Field<NX, NY, NZ, Component::Hz> cze;
};

template <int NX, int NY, int NZ> struct JField {
  using Shape = GridShape<NX, NY, NZ>;
  Field<NX, NY, NZ, Component::Jx> x;
  Field<NX, NY, NZ, Component::Jy> y;
  Field<NX, NY, NZ, Component::Jz> z;
};

template <int NX, int NY, int NZ, Component C>
void setup_field_coeffients(Field<NX, NY, NZ, C> &ca,
                            Field<NX, NY, NZ, C> &cb, double factor) {
  using l = Layout<NX, NY, NZ, C>;
  for (auto x{0}; x < l::nx; ++x)
    for (auto y{0}; y < l::ny; ++y)
      for (auto z{0}; z < l::nz; ++z) {
        ca(x, y, z) = 1.0;
        cb(x, y, z) = factor;
      }
}
