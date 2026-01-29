#include "fdtd_types.hpp"

template <int N> struct DefaultFieldUpdater {
  void half_update_e(EField<N> &E, JField<N> &J, HField<N> &H) const {
    for (auto x{0}; x < E.x.nx(); ++x)
      for (auto y{1}; y < E.x.ny() - 1; ++y)
        for (auto z{1}; z < E.x.nz() - 1; ++z) {
          const double dhz_dy = H.z(x, y, z) - H.z(x, y - 1, z);
          const double dhy_dz = H.y(x, y, z) - H.y(x, y, z - 1);
          const double j = J.x(x, y, z);
          E.x(x, y, z) = E.cxe(x, y, z) * E.x(x, y, z) +
                         0.5 * E.cxh(x, y, z) * (dhz_dy - dhy_dz - j);
        }

    for (auto x{1}; x < E.y.nx() - 1; ++x)
      for (auto y{0}; y < E.y.ny(); ++y)
        for (auto z{1}; z < E.y.nz() - 1; ++z) {
          const double dhx_dhz = H.x(x, y, z) - H.x(x, y, z - 1);
          const double dhz_dhx = H.z(x, y, z) - H.z(x - 1, y, z);
          const double j = J.y(x, y, z);
          E.y(x, y, z) = E.cye(x, y, z) * E.y(x, y, z) +
                         0.5 * E.cyh(x, y, z) * (dhx_dhz - dhz_dhx - j);
        }

    for (auto x{1}; x < E.z.nx() - 1; ++x)
      for (auto y{1}; y < E.z.ny() - 1; ++y)
        for (auto z{0}; z < E.z.nz(); ++z) {
          const double dhy_dhx = H.y(x, y, z) - H.y(x - 1, y, z);
          const double dhx_dhy = H.x(x, y, z) - H.x(x, y - 1, z);
          const double j = J.z(x, y, z);
          E.z(x, y, z) = E.cze(x, y, z) * E.z(x, y, z) +
                         0.5 * E.czh(x, y, z) * (dhy_dhx - dhx_dhy - j);
        }
  }

  void half_update_h(EField<N> &E, HField<N> &H) const {
    for (auto x{0}; x < H.x.nx(); ++x)
      for (auto y{0}; y < H.x.ny(); ++y)
        for (auto z{0}; z < H.x.nz(); ++z) {
          const double dey_dz = E.y(x, y, z + 1) - E.y(x, y, z);
          const double dez_dy = E.z(x, y + 1, z) - E.z(x, y, z);
          H.x(x, y, z) = H.cxh(x, y, z) * H.x(x, y, z) +
                         0.5 * H.cxe(x, y, z) * (dey_dz - dez_dy);
        }

    for (auto x{0}; x < H.y.nx(); ++x)
      for (auto y{0}; y < H.y.ny(); ++y)
        for (auto z{0}; z < H.y.nz(); ++z) {
          const double dex_dz = E.x(x, y, z + 1) - E.x(x, y, z);
          const double dez_dx = E.z(x + 1, y, z) - E.z(x, y, z);
          H.y(x, y, z) = H.cyh(x, y, z) * H.y(x, y, z) +
                         0.5 * H.cye(x, y, z) * (dez_dx - dex_dz);
        }

    for (auto x{0}; x < H.z.nx(); ++x)
      for (auto y{0}; y < H.z.ny(); ++y)
        for (auto z{0}; z < H.z.nz(); ++z) {
          const double dex_dy = E.x(x, y + 1, z) - E.x(x, y, z);
          const double dey_dx = E.y(x + 1, y, z) - E.y(x, y, z);
          H.z(x, y, z) = H.czh(x, y, z) * H.z(x, y, z) +
                         0.5 * H.cze(x, y, z) * (dex_dy - dey_dx);
        }
  }
};

template <int N> struct PMLFieldUpdater {
  int thickness;

  void half_update_e(EField<N> &E, JField<N> &J, HField<N> &H) const {
    for (auto x{0}; x < E.x.nx(); ++x)
      for (auto y{1}; y < E.x.ny() - 1; ++y)
        for (auto z{1}; z < E.x.nz() - 1; ++z) {
          const double dhz_dy = H.z(x, y, z) - H.z(x, y - 1, z);
          const double dhy_dz = H.y(x, y, z) - H.y(x, y, z - 1);
          const double j = J.x(x, y, z);
          E.x(x, y, z) = E.cxe(x, y, z) * E.x(x, y, z) +
                         0.5 * E.cxh(x, y, z) * (dhz_dy - dhy_dz - j);
        }

    for (auto x{1}; x < E.y.nx() - 1; ++x)
      for (auto y{0}; y < E.y.ny(); ++y)
        for (auto z{1}; z < E.y.nz() - 1; ++z) {
          const double dhx_dhz = H.x(x, y, z) - H.x(x, y, z - 1);
          const double dhz_dhx = H.z(x, y, z) - H.z(x - 1, y, z);
          const double j = J.y(x, y, z);
          E.y(x, y, z) = E.cye(x, y, z) * E.y(x, y, z) +
                         0.5 * E.cyh(x, y, z) * (dhx_dhz - dhz_dhx - j);
        }

    for (auto x{1}; x < E.z.nx() - 1; ++x)
      for (auto y{1}; y < E.z.ny() - 1; ++y)
        for (auto z{0}; z < E.z.nz(); ++z) {
          const double dhy_dhx = H.y(x, y, z) - H.y(x - 1, y, z);
          const double dhx_dhy = H.x(x, y, z) - H.x(x, y - 1, z);
          const double j = J.z(x, y, z);
          E.z(x, y, z) = E.cze(x, y, z) * E.z(x, y, z) +
                         0.5 * E.czh(x, y, z) * (dhy_dhx - dhx_dhy - j);
        }
  }

  void half_update_h(EField<N> &E, HField<N> &H) const {
    for (auto x{0}; x < H.x.nx(); ++x)
      for (auto y{0}; y < H.x.ny(); ++y)
        for (auto z{0}; z < H.x.nz(); ++z) {
          const double dey_dz = E.y(x, y, z + 1) - E.y(x, y, z);
          const double dez_dy = E.z(x, y + 1, z) - E.z(x, y, z);
          H.x(x, y, z) = H.cxh(x, y, z) * H.x(x, y, z) +
                         0.5 * H.cxe(x, y, z) * (dey_dz - dez_dy);
        }

    for (auto x{0}; x < H.y.nx(); ++x)
      for (auto y{0}; y < H.y.ny(); ++y)
        for (auto z{0}; z < H.y.nz(); ++z) {
          const double dex_dz = E.x(x, y, z + 1) - E.x(x, y, z);
          const double dez_dx = E.z(x + 1, y, z) - E.z(x, y, z);
          H.y(x, y, z) = H.cyh(x, y, z) * H.y(x, y, z) +
                         0.5 * H.cye(x, y, z) * (dez_dx - dex_dz);
        }

    for (auto x{0}; x < H.z.nx(); ++x)
      for (auto y{0}; y < H.z.ny(); ++y)
        for (auto z{0}; z < H.z.nz(); ++z) {
          const double dex_dy = E.x(x, y + 1, z) - E.x(x, y, z);
          const double dey_dx = E.y(x + 1, y, z) - E.y(x, y, z);
          H.z(x, y, z) = H.czh(x, y, z) * H.z(x, y, z) +
                         0.5 * H.cze(x, y, z) * (dex_dy - dey_dx);
        }
  }
};
