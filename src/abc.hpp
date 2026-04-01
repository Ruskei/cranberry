#pragma once

#include <vector>

#include "config.hpp"
#include "fdtd_types.hpp"

template <int NX, int NY, int NZ> struct ABC {
  using Shape = GridShape<NX, NY, NZ>;
  template <int SX, int SY> struct Slab {
    std::vector<double> v;

    Slab() : v(SX * SY) {}

    double &operator()(int a, int b) { return v[a * SY + b]; }
    const double &operator()(int a, int b) const { return v[a * SY + b]; }
  };

  using XWallEy = Slab<NY - 1, NZ>;
  using XWallEz = Slab<NY, NZ - 1>;
  using YWallEx = Slab<NX - 1, NZ>;
  using YWallEz = Slab<NX, NZ - 1>;
  using ZWallEx = Slab<NX - 1, NY>;
  using ZWallEy = Slab<NX, NY - 1>;

  XWallEy eyx0;
  XWallEz ezx0;
  XWallEy eyx1;
  XWallEz ezx1;

  YWallEx exy0;
  YWallEz ezy0;
  YWallEx exy1;
  YWallEz ezy1;

  ZWallEx exz0;
  ZWallEy eyz0;
  ZWallEx exz1;
  ZWallEy eyz1;

  void apply(EField<NX, NY, NZ> &E) {
    const double abcco = Config::abcco;
    auto &ex = E.x;
    auto &ey = E.y;
    auto &ez = E.z;

    {
      const int x{0};
      for (int y{0}; y < NY - 1; ++y)
        for (int z{0}; z < NZ; ++z) {
          ey(x, y, z) = eyx0(y, z) + abcco * (ey(x + 1, y, z) - ey(x, y, z));
          eyx0(y, z) = ey(x + 1, y, z);
        }

      for (int y{0}; y < NY; ++y)
        for (int z{0}; z < NZ - 1; ++z) {
          const double co = ezx0(y, z);
          ez(x, y, z) = co + abcco * (ez(x + 1, y, z) - ez(x, y, z));
          ezx0(y, z) = ez(x + 1, y, z);
        }
    }

    {
      const int x = NX - 1;
      for (int y{0}; y < NY - 1; ++y)
        for (int z{0}; z < NZ; ++z) {
          ey(x, y, z) = eyx1(y, z) + abcco * (ey(x - 1, y, z) - ey(x, y, z));
          eyx1(y, z) = ey(x - 1, y, z);
        }

      for (int y{0}; y < NY; ++y)
        for (int z{0}; z < NZ - 1; ++z) {
          ez(x, y, z) = ezx1(y, z) + abcco * (ez(x - 1, y, z) - ez(x, y, z));
          ezx1(y, z) = ez(x - 1, y, z);
        }
    }

    {
      const int y{0};
      for (int x{0}; x < NX - 1; ++x)
        for (int z{0}; z < NZ; ++z) {
          ex(x, y, z) = exy0(x, z) + abcco * (ex(x, y + 1, z) - ex(x, y, z));
          exy0(x, z) = ex(x, y + 1, z);
        }

      for (int x{0}; x < NX; ++x)
        for (int z{0}; z < NZ - 1; ++z) {
          ez(x, y, z) = ezy0(x, z) + abcco * (ez(x, y + 1, z) - ez(x, y, z));
          ezy0(x, z) = ez(x, y + 1, z);
        }
    }

    {
      const int y = NY - 1;
      for (int x{0}; x < NX - 1; ++x)
        for (int z{0}; z < NZ; ++z) {
          ex(x, y, z) = exy1(x, z) + abcco * (ex(x, y - 1, z) - ex(x, y, z));
          exy1(x, z) = ex(x, y - 1, z);
        }

      for (int x{0}; x < NX; ++x)
        for (int z{0}; z < NZ - 1; ++z) {
          ez(x, y, z) = ezy1(x, z) + abcco * (ez(x, y - 1, z) - ez(x, y, z));
          ezy1(x, z) = ez(x, y - 1, z);
        }
    }

    {
      const int z{0};
      for (int x{0}; x < NX - 1; ++x)
        for (int y{0}; y < NY; ++y) {
          ex(x, y, z) = exz0(x, y) + abcco * (ex(x, y, z + 1) - ex(x, y, z));
          exz0(x, y) = ex(x, y, z + 1);
        }

      for (int x{0}; x < NX; ++x)
        for (int y{0}; y < NY - 1; ++y) {
          ey(x, y, z) = eyz0(x, y) + abcco * (ey(x, y, z + 1) - ey(x, y, z));
          eyz0(x, y) = ey(x, y, z + 1);
        }
    }

    {
      const int z = NZ - 1;
      for (int x{0}; x < NX - 1; ++x)
        for (int y{0}; y < NY; ++y) {
          ex(x, y, z) = exz1(x, y) + abcco * (ex(x, y, z - 1) - ex(x, y, z));
          exz1(x, y) = ex(x, y, z - 1);
        }

      for (int x{0}; x < NX; ++x)
        for (int y{0}; y < NY - 1; ++y) {
          ey(x, y, z) = eyz1(x, y) + abcco * (ey(x, y, z - 1) - ey(x, y, z));
          eyz1(x, y) = ey(x, y, z - 1);
        }
    }
  }
};
