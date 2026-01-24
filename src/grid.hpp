#pragma once

#include "config.hpp"
#include "fdtd_types.hpp"

template <int N> struct Grid {
  Field<N, Component::Ex> ex;
  Field<N, Component::Ex> cexe;
  Field<N, Component::Ex> cexh;

  Field<N, Component::Ey> ey;
  Field<N, Component::Ey> ceye;
  Field<N, Component::Ey> ceyh;

  Field<N, Component::Ez> ez;
  Field<N, Component::Ez> ceze;
  Field<N, Component::Ez> cezh;

  Field<N, Component::Jx> jx;
  Field<N, Component::Jy> jy;
  Field<N, Component::Jz> jz;

  Field<N, Component::Hx> hx;
  Field<N, Component::Hx> chxh;
  Field<N, Component::Hx> chxe;

  Field<N, Component::Hy> hy;
  Field<N, Component::Hy> chyh;
  Field<N, Component::Hy> chye;

  Field<N, Component::Hz> hz;
  Field<N, Component::Hz> chzh;
  Field<N, Component::Hz> chze;

  void setup_coefficients() {
    const double c = Config::courants;
    const double half_c = c * 0.5; // simulate in half-steps
    const double imp0 = Config::imp0;
    setup_field_coeffients(cexe, cexh, half_c * imp0);
    setup_field_coeffients(ceye, ceyh, half_c * imp0);
    setup_field_coeffients(ceze, cezh, half_c * imp0);
    setup_field_coeffients(chxh, chxe, half_c / imp0);
    setup_field_coeffients(chyh, chye, half_c / imp0);
    setup_field_coeffients(chzh, chze, half_c / imp0);
  }

  void update_h() {
    for (auto x{0}; x < hx.nx(); ++x)
      for (auto y{0}; y < hx.ny(); ++y)
        for (auto z{0}; z < hx.nz(); ++z) {
          const double dey_dz = ey(x, y, z + 1) - ey(x, y, z);
          const double dez_dy = ez(x, y + 1, z) - ez(x, y, z);
          hx(x, y, z) =
              chxh(x, y, z) * hx(x, y, z) + chxe(x, y, z) * (dey_dz - dez_dy);
        }

    for (auto x{0}; x < hy.nx(); ++x)
      for (auto y{0}; y < hy.ny(); ++y)
        for (auto z{0}; z < hy.nz(); ++z) {
          const double dex_dz = ex(x, y, z + 1) - ex(x, y, z);
          const double dez_dx = ez(x + 1, y, z) - ez(x, y, z);
          hy(x, y, z) =
              chyh(x, y, z) * hy(x, y, z) + chye(x, y, z) * (dez_dx - dex_dz);
        }

    for (auto x{0}; x < hz.nx(); ++x)
      for (auto y{0}; y < hz.ny(); ++y)
        for (auto z{0}; z < hz.nz(); ++z) {
          const double dex_dy = ex(x, y + 1, z) - ex(x, y, z);
          const double dey_dx = ey(x + 1, y, z) - ey(x, y, z);
          hz(x, y, z) =
              chzh(x, y, z) * hz(x, y, z) + chze(x, y, z) * (dex_dy - dey_dx);
        }
  }

  void update_e() {
    for (auto x{0}; x < ex.nx(); ++x)
      for (auto y{1}; y < ex.ny() - 1; ++y)
        for (auto z{1}; z < ex.nz() - 1; ++z) {
          const double dhz_dy = hz(x, y, z) - hz(x, y - 1, z);
          const double dhy_dz = hy(x, y, z) - hy(x, y, z - 1);
          const double j = jx(x, y, z);
          ex(x, y, z) = cexe(x, y, z) * ex(x, y, z) +
                        cexh(x, y, z) * (dhz_dy - dhy_dz - j);
        }

    for (auto x{1}; x < ey.nx() - 1; ++x)
      for (auto y{0}; y < ey.ny(); ++y)
        for (auto z{1}; z < ey.nz() - 1; ++z) {
          const double dhx_dhz = hx(x, y, z) - hx(x, y, z - 1);
          const double dhz_dhx = hz(x, y, z) - hz(x - 1, y, z);
          const double j = jy(x, y, z);
          ey(x, y, z) = ceye(x, y, z) * ey(x, y, z) +
                        ceyh(x, y, z) * (dhx_dhz - dhz_dhx - j);
        }

    for (auto x{1}; x < ez.nx() - 1; ++x)
      for (auto y{1}; y < ez.ny() - 1; ++y)
        for (auto z{0}; z < ez.nz(); ++z) {
          const double dhy_dhx = hy(x, y, z) - hy(x - 1, y, z);
          const double dhx_dhy = hx(x, y, z) - hx(x, y - 1, z);
          const double j = jz(x, y, z);
          ez(x, y, z) = ceze(x, y, z) * ez(x, y, z) +
                        cezh(x, y, z) * (dhy_dhx - dhx_dhy - j);
        }
  }
};
