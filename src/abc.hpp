#pragma once

#include <vector>

#include "grid.hpp"

template <int N> struct ABC {
  static constexpr int size = (N - 1) * N;

  template <int NX, int NY> struct Slab {
    std::vector<double> v;

    Slab() : v(NX * NY) {}

    double &operator()(int a, int b) { return v[a * NY + b]; }
    const double &operator()(int a, int b) const { return v[a * NY + b]; }
  };

  using WallA = Slab<N - 1, N>;
  using WallB = Slab<N, N - 1>;

  WallA eyx0;
  WallB ezx0;
  WallA eyx1;
  WallB ezx1;

  WallA exy0;
  WallB ezy0;
  WallA exy1;
  WallB ezy1;

  WallA exz0;
  WallB eyz0;
  WallA exz1;
  WallB eyz1;

  void apply(Grid<N> &grid) {
    const double abcco = Config::abcco;
    auto &ex = grid.ex;
    auto &ey = grid.ey;
    auto &ez = grid.ez;

    {
      const int x{0};
      for (int y{0}; y < N - 1; ++y)
        for (int z{0}; z < N; ++z) {
          ey(x, y, z) = eyx0(y, z) + abcco * (ey(x + 1, y, z) - ey(x, y, z));
          eyx0(y, z) = ey(x + 1, y, z);
        }

      for (int y{0}; y < N; ++y)
        for (int z{0}; z < N - 1; ++z) {
          const double co = ezx0(y, z);
          ez(x, y, z) = co + abcco * (ez(x + 1, y, z) - ez(x, y, z));
          ezx0(y, z) = ez(x + 1, y, z);
        }
    }

    {
      const int x = N - 1;
      for (int y{0}; y < N - 1; ++y)
        for (int z{0}; z < N; ++z) {
          ey(x, y, z) = eyx1(y, z) + abcco * (ey(x - 1, y, z) - ey(x, y, z));
          eyx1(y, z) = ey(x - 1, y, z);
        }

      for (int y{0}; y < N; ++y)
        for (int z{0}; z < N - 1; ++z) {
          ez(x, y, z) = ezx1(y, z) + abcco * (ez(x - 1, y, z) - ez(x, y, z));
          ezx1(y, z) = ez(x - 1, y, z);
        }
    }

    {
      const int y{0};
      for (int x{0}; x < N - 1; ++x)
        for (int z{0}; z < N; ++z) {
          ex(x, y, z) = exy0(x, z) + abcco * (ex(x, y + 1, z) - ex(x, y, z));
          exy0(x, z) = ex(x, y + 1, z);
        }

      for (int x{0}; x < N; ++x)
        for (int z{0}; z < N - 1; ++z) {
          ez(x, y, z) = ezy0(x, z) + abcco * (ez(x, y + 1, z) - ez(x, y, z));
          ezy0(x, z) = ez(x, y + 1, z);
        }
    }

    {
      const int y = N - 1;
      for (int x{0}; x < N - 1; ++x)
        for (int z{0}; z < N; ++z) {
          ex(x, y, z) = exy1(x, z) + abcco * (ex(x, y - 1, z) - ex(x, y, z));
          exy1(x, z) = ex(x, y - 1, z);
        }

      for (int x{0}; x < N; ++x)
        for (int z{0}; z < N - 1; ++z) {
          ez(x, y, z) = ezy1(x, z) + abcco * (ez(x, y - 1, z) - ez(x, y, z));
          ezy1(x, z) = ez(x, y - 1, z);
        }
    }

    {
      const int z{0};
      for (int x{0}; x < N - 1; ++x)
        for (int y{0}; y < N; ++y) {
          ex(x, y, z) = exz0(x, y) + abcco * (ex(x, y, z + 1) - ex(x, y, z));
          exz0(x, y) = ex(x, y, z + 1);
        }

      for (int x{0}; x < N; ++x)
        for (int y{0}; y < N - 1; ++y) {
          ey(x, y, z) = eyz0(x, y) + abcco * (ey(x, y, z + 1) - ey(x, y, z));
          eyz0(x, y) = ey(x, y, z + 1);
        }
    }

    {
      const int z = N - 1;
      for (int x{0}; x < N - 1; ++x)
        for (int y{0}; y < N; ++y) {
          ex(x, y, z) = exz1(x, y) + abcco * (ex(x, y, z - 1) - ex(x, y, z));
          exz1(x, y) = ex(x, y, z - 1);
        }

      for (int x{0}; x < N; ++x)
        for (int y{0}; y < N - 1; ++y) {
          ey(x, y, z) = eyz1(x, y) + abcco * (ey(x, y, z - 1) - ey(x, y, z));
          eyz1(x, y) = ey(x, y, z - 1);
        }
    }
  }
};
