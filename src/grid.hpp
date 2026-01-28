#pragma once

#include <algorithm>
#include <iostream>
#include <utility>

#include "config.hpp"
#include "electrostatics.hpp"
#include "fdtd_types.hpp"
#include "particles.hpp"

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

  Field<N, Component::Charge> charge;
  Field<N, Component::Potential> potential;

  std::vector<Particle> particles;

  void deposit_charge(const std::vector<Particle> &particles,
                      Field<N, Component::Charge> &charge) {
    for (const auto &p : particles) {
      const double px = p.px;
      const double py = p.py;
      const double pz = p.pz;

      const int i = static_cast<int>(std::floor(px));
      const int j = static_cast<int>(std::floor(py));
      const int k = static_cast<int>(std::floor(pz));

      const int is = std::max(0, i - 1);
      const int js = std::max(0, j - 1);
      const int ks = std::max(0, k - 1);

      const int ie = std::min(charge.nx() - 1, i + 1);
      const int je = std::min(charge.ny() - 1, j + 1);
      const int ke = std::min(charge.nz() - 1, k + 1);

      for (auto x{is}; x <= ie; ++x)
        for (auto y{js}; y <= je; ++y)
          for (auto z{ks}; z <= ke; ++z) {
            charge(x, y, z) += -p.q * form_factor(px, x) * form_factor(py, y) *
                               form_factor(pz, z);
          }
    }
  }

  void apply_potential() {
    for (auto x{1}; x < ex.nx() - 1; ++x)
      for (auto y{2}; y < ex.ny() - 2; ++y)
        for (auto z{2}; z < ex.nz() - 2; ++z)
          ex(x, y, z) = potential(x + 1, y, z) - potential(x, y, z);

    for (auto x{2}; x < ey.nx() - 2; ++x)
      for (auto y{1}; y < ey.ny() - 1; ++y)
        for (auto z{2}; z < ey.nz() - 2; ++z)
          ey(x, y, z) = potential(x, y + 1, z) - potential(x, y, z);

    for (auto x{2}; x < ez.nx() - 2; ++x)
      for (auto y{2}; y < ez.ny() - 2; ++y)
        for (auto z{1}; z < ez.nz() - 1; ++z)
          ez(x, y, z) = potential(x, y, z + 1) - potential(x, y, z);
  }

  Grid(std::vector<Particle> particles) : particles{std::move(particles)} {
    std::cout << "Depositing charge..." << std::endl;
    deposit_charge(this->particles, charge);
    binom_smooth_field(charge);
    std::cout << "Solving potential" << std::endl;
    solve_potential(charge, potential);
    std::cout << "Applying potential" << std::endl;
    apply_potential();
    std::cout << "Finished Grid initialization" << std::endl;
  }

  void setup_coefficients() {
    const double c = Config::courants;
    const double imp0 = Config::imp0;
    setup_field_coeffients(cexe, cexh, c * imp0);
    setup_field_coeffients(ceye, ceyh, c * imp0);
    setup_field_coeffients(ceze, cezh, c * imp0);
    setup_field_coeffients(chxh, chxe, c / imp0);
    setup_field_coeffients(chyh, chye, c / imp0);
    setup_field_coeffients(chzh, chze, c / imp0);
  }

  void half_update_h() {
    for (auto x{0}; x < hx.nx(); ++x)
      for (auto y{0}; y < hx.ny(); ++y)
        for (auto z{0}; z < hx.nz(); ++z) {
          const double dey_dz = ey(x, y, z + 1) - ey(x, y, z);
          const double dez_dy = ez(x, y + 1, z) - ez(x, y, z);
          hx(x, y, z) = chxh(x, y, z) * hx(x, y, z) +
                        0.5 * chxe(x, y, z) * (dey_dz - dez_dy);
        }

    for (auto x{0}; x < hy.nx(); ++x)
      for (auto y{0}; y < hy.ny(); ++y)
        for (auto z{0}; z < hy.nz(); ++z) {
          const double dex_dz = ex(x, y, z + 1) - ex(x, y, z);
          const double dez_dx = ez(x + 1, y, z) - ez(x, y, z);
          hy(x, y, z) = chyh(x, y, z) * hy(x, y, z) +
                        0.5 * chye(x, y, z) * (dez_dx - dex_dz);
        }

    for (auto x{0}; x < hz.nx(); ++x)
      for (auto y{0}; y < hz.ny(); ++y)
        for (auto z{0}; z < hz.nz(); ++z) {
          const double dex_dy = ex(x, y + 1, z) - ex(x, y, z);
          const double dey_dx = ey(x + 1, y, z) - ey(x, y, z);
          hz(x, y, z) = chzh(x, y, z) * hz(x, y, z) +
                        0.5 * chze(x, y, z) * (dex_dy - dey_dx);
        }
  }

  void half_update_e() {
    for (auto x{0}; x < ex.nx(); ++x)
      for (auto y{1}; y < ex.ny() - 1; ++y)
        for (auto z{1}; z < ex.nz() - 1; ++z) {
          const double dhz_dy = hz(x, y, z) - hz(x, y - 1, z);
          const double dhy_dz = hy(x, y, z) - hy(x, y, z - 1);
          const double j = jx(x, y, z);
          ex(x, y, z) = cexe(x, y, z) * ex(x, y, z) +
                        0.5 * cexh(x, y, z) * (dhz_dy - dhy_dz - j);
        }

    for (auto x{1}; x < ey.nx() - 1; ++x)
      for (auto y{0}; y < ey.ny(); ++y)
        for (auto z{1}; z < ey.nz() - 1; ++z) {
          const double dhx_dhz = hx(x, y, z) - hx(x, y, z - 1);
          const double dhz_dhx = hz(x, y, z) - hz(x - 1, y, z);
          const double j = jy(x, y, z);
          ey(x, y, z) = ceye(x, y, z) * ey(x, y, z) +
                        0.5 * ceyh(x, y, z) * (dhx_dhz - dhz_dhx - j);
        }

    for (auto x{1}; x < ez.nx() - 1; ++x)
      for (auto y{1}; y < ez.ny() - 1; ++y)
        for (auto z{0}; z < ez.nz(); ++z) {
          const double dhy_dhx = hy(x, y, z) - hy(x - 1, y, z);
          const double dhx_dhy = hx(x, y, z) - hx(x, y - 1, z);
          const double j = jz(x, y, z);
          ez(x, y, z) = ceze(x, y, z) * ez(x, y, z) +
                        0.5 * cezh(x, y, z) * (dhy_dhx - dhx_dhy - j);
        }
  }

  void step_particles() {
    for (auto &p : particles) p.age++;
  }

  void push_particles() {
    const double epsilon = 1e-12;
    const double dt = Config::dt;
    for (auto &p : particles) {
      const double nx = p.px;
      const double ny = p.py;
      const double nz = p.pz;

      double vx = p.vx;
      double vy = p.vy;
      double vz = p.vz;

      double vv = vx * vx + vy * vy + vz * vz;
      if (!std::isfinite(vv) || vv >= 1.0 - epsilon) {
        double scale = std::sqrt((1.0 - 1e-12) / std::max(vv, epsilon));
        vx *= scale;
        vy *= scale;
        vz *= scale;
        vv = vx * vx + vy * vy + vz * vz;
      }

      const double q = p.q;
      const double m = p.m;

      double lorentz = 1.0 / std::sqrt(1.0 - vv);

      const double exv = ex.interpolate_at(nx, ny, nz);
      const double eyv = ey.interpolate_at(nx, ny, nz);
      const double ezv = ez.interpolate_at(nx, ny, nz);

      const double hxv = hx.interpolate_at(nx, ny, nz);
      const double hyv = hy.interpolate_at(nx, ny, nz);
      const double hzv = hz.interpolate_at(nx, ny, nz);

      const double ux = lorentz * vx;
      const double uy = lorentz * vy;
      const double uz = lorentz * vz;

      const double u_minus_x = ux + (q * dt / (2.0 * m)) * exv;
      const double u_minus_y = uy + (q * dt / (2.0 * m)) * eyv;
      const double u_minus_z = uz + (q * dt / (2.0 * m)) * ezv;
      const double umum = (u_minus_x * u_minus_x + u_minus_y * u_minus_y +
                           u_minus_z * u_minus_z);

      const double lorentz_bar = 1.0 / sqrt(1.0 + umum);
      const double tx = q * dt / (2.0 * m) * hxv * lorentz_bar;
      const double ty = q * dt / (2.0 * m) * hyv * lorentz_bar;
      const double tz = q * dt / (2.0 * m) * hzv * lorentz_bar;
      const double tt = tx * tx + ty * ty + tz * tz;

      const double u_prime_x = u_minus_x + (u_minus_y * tz - u_minus_z * ty);
      const double u_prime_y = u_minus_y + (u_minus_z * tx - u_minus_x * tz);
      const double u_prime_z = u_minus_z + (u_minus_x * ty - u_minus_y * tx);

      const double ax = 2 * tx / (1 + tt);
      const double ay = 2 * ty / (1 + tt);
      const double az = 2 * tz / (1 + tt);

      const double u_plus_x = u_minus_x + (u_prime_y * az - u_prime_z * ay);
      const double u_plus_y = u_minus_y + (u_prime_z * ax - u_prime_x * az);
      const double u_plus_z = u_minus_z + (u_prime_x * ay - u_prime_y * ax);

      const double u_next_x = u_plus_x + q * dt / (2 * m) * exv;
      const double u_next_y = u_plus_y + q * dt / (2 * m) * eyv;
      const double u_next_z = u_plus_z + q * dt / (2 * m) * ezv;
      const double unun =
          (u_next_x * u_next_x + u_next_y * u_next_y + u_next_z * u_next_z);

      const double lorentz_next = std::sqrt(1.0 + unun);
      p.vx = u_next_x / lorentz_next;
      p.vy = u_next_y / lorentz_next;
      p.vz = u_next_z / lorentz_next;

      p.p_prev_x = p.px;
      p.p_prev_y = p.py;
      p.p_prev_z = p.pz;

      p.px += p.vx * dt;
      p.py += p.vy * dt;
      p.pz += p.vz * dt;
    }
  }

  void deposit_currents() {
    std::fill(jx.v.begin(), jx.v.end(), 0);
    std::fill(jy.v.begin(), jy.v.end(), 0);
    std::fill(jz.v.begin(), jz.v.end(), 0);

    for (const auto &p : particles) {
      const double px = p.px;
      const double py = p.py;
      const double pz = p.pz;

      const double ppx = p.p_prev_x;
      const double ppy = p.p_prev_y;
      const double ppz = p.p_prev_z;

      const int i = static_cast<int>(std::floor(px));
      const int j = static_cast<int>(std::floor(py));
      const int k = static_cast<int>(std::floor(pz));

      const int pi = static_cast<int>(std::floor(ppx));
      const int pj = static_cast<int>(std::floor(ppy));
      const int pk = static_cast<int>(std::floor(ppz));

      const int is = std::max(0, std::min(i, pi) - 1);
      const int js = std::max(0, std::min(j, pj) - 1);
      const int ks = std::max(0, std::min(k, pk) - 1);

      const int ie = std::min(jx.nx() - 1, std::max(i, pi) + 1);
      const int je = std::min(jy.ny() - 1, std::max(j, pj) + 1);
      const int ke = std::min(jz.nz() - 1, std::max(k, pk) + 1);

      const double move_co = -p.q / Config::dt;

      for (auto y{js}; y <= je; ++y)
        for (auto z{ks}; z <= ke; ++z) {
          double move_sum = 0.0;
          for (auto x{is}; x <= ie; ++x) {
            move_sum +=
                move_co * form_factor_diff<CartesianComponent::x>(p, x, y, z);
            jx(x, y, z) += move_sum;
          }
        }

      for (auto x{is}; x <= ie; ++x)
        for (auto z{ks}; z <= ke; ++z) {
          double move_sum = 0.0;
          for (auto y{js}; y <= je; ++y) {
            move_sum +=
                move_co * form_factor_diff<CartesianComponent::y>(p, x, y, z);
            jy(x, y, z) += move_sum;
          }
        }

      for (auto x{is}; x <= ie; ++x)
        for (auto y{js}; y <= je; ++y) {
          double move_sum = 0.0;
          for (auto z{ks}; z <= ke; ++z) {
            move_sum +=
                move_co * form_factor_diff<CartesianComponent::z>(p, x, y, z);
            jz(x, y, z) += move_sum;
          }
        }
    }
  }

  template<Component C>
  void binom_smooth_field(Field<N, C> &field) {
    for (auto y{0}; y < field.ny(); ++y)
      for (auto z{0}; z < field.nz(); ++z) {
        double left = field(0, y, z);
        double center = field(1, y, z);
        for (auto x{1}; x < field.nx() - 1; ++x) {
          const double right = field(x + 1, y, z);
          field(x, y, z) = 0.25 * left + 0.5 * center + 0.25 * right;

          left = center;
          center = right;
        }
      }
    for (auto x{0}; x < field.nx(); ++x)
      for (auto z{0}; z < field.nz(); ++z) {
        double left = field(x, 0, z);
        double center = field(x, 1, z);
        for (auto y{1}; y < field.ny() - 1; ++y) {
          const double right = field(x, y + 1, z);
          field(x, y, z) = 0.25 * left + 0.5 * center + 0.25 * right;

          left = center;
          center = right;
        }
      }
    for (auto x{0}; x < field.nx(); ++x)
      for (auto y{0}; y < field.ny(); ++y) {
        double left = field(x, y, 0);
        double center = field(x, y, 1);
        for (auto z{1}; z < field.nz() - 1; ++z) {
          const double right = field(x, y, z + 1);
          field(x, y, z) = 0.25 * left + 0.5 * center + 0.25 * right;

          left = center;
          center = right;
        }
      }
  }

  void smooth_currents() {
    binom_smooth_field(jx);
    binom_smooth_field(jy);
    binom_smooth_field(jz);
  }

  void check_currents() {
    for (const auto &p : particles) {
      double residual{};

      const double px = p.px;
      const double py = p.py;
      const double pz = p.pz;

      const double ppx = p.p_prev_x;
      const double ppy = p.p_prev_y;
      const double ppz = p.p_prev_z;

      const int i = static_cast<int>(std::floor(px));
      const int j = static_cast<int>(std::floor(py));
      const int k = static_cast<int>(std::floor(pz));

      const int pi = static_cast<int>(std::floor(ppx));
      const int pj = static_cast<int>(std::floor(ppy));
      const int pk = static_cast<int>(std::floor(ppz));

      const int is = std::max(0, std::min(i, pi) - 1);
      const int js = std::max(0, std::min(j, pj) - 1);
      const int ks = std::max(0, std::min(k, pk) - 1);

      const int ie = std::min(jx.nx() - 1, std::max(i, pi) + 1);
      const int je = std::min(jy.ny() - 1, std::max(j, pj) + 1);
      const int ke = std::min(jz.nz() - 1, std::max(k, pk) + 1);

      for (auto x{is}; x <= ie; ++x)
        for (auto y{js}; y <= je; ++y)
          for (auto z{ks}; z <= ke; ++z) {
            const double charge = p.q * form_factor(px, x) *
                                  form_factor(py, y) * form_factor(pz, z);
            const double charge_prev = p.q * form_factor(ppx, x) *
                                       form_factor(ppy, y) *
                                       form_factor(ppz, z);

            const double clamped_x = std::max(0, x - 1);
            const double clamped_y = std::max(0, y - 1);
            const double clamped_z = std::max(0, z - 1);

            const double div_J = (jx(x, y, z) - jx(clamped_x, y, z)) +
                                 (jy(x, y, z) - jy(x, clamped_y, z)) +
                                 (jz(x, y, z) - jz(x, y, clamped_z));

            // if (std::abs(div_J) > 0.0)
            //   std::cout << "∇∙J_(" << x << "," << y << "," << z << ")=" <<
            //   div_J << std::endl;

            residual += ((charge - charge_prev) / Config::dt + div_J);
          }

      if (std::abs(residual) > 1e-11)
        std::cout << "continuity residual: " << residual << " at age: " << p.age
                  << std::endl;
    }
  }
};
