#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <utility>
#include <vector>

#include "config.hpp"
#include "fdtd_types.hpp"
#include "particles.hpp"

template <int NX, int NY, int NZ> struct GaussLawResult {
  static constexpr std::size_t metric_count{8};

  std::string name;
  int interval{1};
  std::vector<std::array<double, metric_count>> metrics;

  GaussLawResult(std::string name, int interval)
      : name{std::move(name)}, interval{interval} {
    mkdir("out", 0755);
    mkdir(("out/" + this->name).c_str(), 0755);
  }

  template <Component C>
  static void binom_smooth_field(Field<NX, NY, NZ, C> &field) {
    for (int y = 0; y < field.ny(); ++y)
      for (int z = 0; z < field.nz(); ++z) {
        double left = field(0, y, z);
        double center = field(1, y, z);
        for (int x = 1; x < field.nx() - 1; ++x) {
          const double right = field(x + 1, y, z);
          field(x, y, z) = 0.25 * left + 0.5 * center + 0.25 * right;
          left = center;
          center = right;
        }
      }

    for (int x = 0; x < field.nx(); ++x)
      for (int z = 0; z < field.nz(); ++z) {
        double left = field(x, 0, z);
        double center = field(x, 1, z);
        for (int y = 1; y < field.ny() - 1; ++y) {
          const double right = field(x, y + 1, z);
          field(x, y, z) = 0.25 * left + 0.5 * center + 0.25 * right;
          left = center;
          center = right;
        }
      }

    for (int x = 0; x < field.nx(); ++x)
      for (int y = 0; y < field.ny(); ++y) {
        double left = field(x, y, 0);
        double center = field(x, y, 1);
        for (int z = 1; z < field.nz() - 1; ++z) {
          const double right = field(x, y, z + 1);
          field(x, y, z) = 0.25 * left + 0.5 * center + 0.25 * right;
          left = center;
          center = right;
        }
      }
  }

  void step(double time, const EField<NX, NY, NZ> &E,
            const std::vector<Particle> &particles) {
    if (static_cast<int>(std::round(time / Config::dt)) % interval != 0)
      return;

    Field<NX, NY, NZ, Component::Charge> live_charge;
    for (const auto &p : particles) {
      const int i = static_cast<int>(std::floor(p.px));
      const int j = static_cast<int>(std::floor(p.py));
      const int k = static_cast<int>(std::floor(p.pz));

      const int is = std::max(0, i - particle_radius);
      const int js = std::max(0, j - particle_radius);
      const int ks = std::max(0, k - particle_radius);
      const int ie = std::min(live_charge.nx() - 1, i + 1 + particle_radius);
      const int je = std::min(live_charge.ny() - 1, j + 1 + particle_radius);
      const int ke = std::min(live_charge.nz() - 1, k + 1 + particle_radius);

      if (is > ie || js > je || ks > ke)
        continue;

      for (int x = is; x <= ie; ++x)
        for (int y = js; y <= je; ++y)
          for (int z = ks; z <= ke; ++z)
            live_charge(x, y, z) +=
                p.q * form_factor(p.px, x) * form_factor(p.py, y) *
                form_factor(p.pz, z);
    }

    binom_smooth_field(live_charge);

    double total_div_E = 0.0;
    double total_charge = 0.0;
    double residual_l2_sq = 0.0;
    double charge_l2_sq = 0.0;
    double residual_linf = 0.0;
    double charge_linf = 0.0;

    for (int x = 2; x < live_charge.nx() - 2; ++x)
      for (int y = 2; y < live_charge.ny() - 2; ++y)
        for (int z = 2; z < live_charge.nz() - 2; ++z) {
          const double div_E =
              (E.x(x, y, z) - E.x(x - 1, y, z)) +
              (E.y(x, y, z) - E.y(x, y - 1, z)) +
              (E.z(x, y, z) - E.z(x, y, z - 1));
          const double rho = live_charge(x, y, z);
          const double residual = div_E - rho;

          total_div_E += div_E;
          total_charge += rho;
          residual_l2_sq += residual * residual;
          charge_l2_sq += rho * rho;
          residual_linf = std::max(residual_linf, std::abs(residual));
          charge_linf = std::max(charge_linf, std::abs(rho));
        }

    const double residual_l2 = std::sqrt(residual_l2_sq);
    const double charge_l2 = std::sqrt(charge_l2_sq);
    const double undefined = std::numeric_limits<double>::quiet_NaN();
    const double relative_l2 =
        charge_l2 > 0.0 ? 100.0 * residual_l2 / charge_l2 : undefined;
    const double relative_linf =
        charge_linf > 0.0 ? 100.0 * residual_linf / charge_linf : undefined;

    metrics.push_back({time, total_div_E, total_charge,
                       total_div_E - total_charge, residual_l2, relative_l2,
                       residual_linf, relative_linf});
  }

  void finish() const {
    std::ofstream columns("out/" + name + "/columns.txt");
    columns << "time\n"
               "sum_div_E\n"
               "sum_rho\n"
               "sum_residual\n"
               "residual_L2\n"
               "relative_L2_percent\n"
               "residual_Linf\n"
               "relative_Linf_percent\n";

    std::ofstream output("out/" + name + "/metrics.npy", std::ios::binary);
    const std::string dictionary =
        "{'descr': '<f8', 'fortran_order': False, 'shape': (" +
        std::to_string(metrics.size()) + ", " +
        std::to_string(metric_count) + "), }";
    std::string header = dictionary;
    const std::size_t padding = 64 - ((10 + header.size() + 1) % 64);
    header.append(padding, ' ');
    header.push_back('\n');

    output.write("\x93NUMPY", 6);
    output.put(1);
    output.put(0);
    const auto header_size = static_cast<std::uint16_t>(header.size());
    output.put(static_cast<char>(header_size & 0xff));
    output.put(static_cast<char>((header_size >> 8) & 0xff));
    output.write(header.data(), static_cast<std::streamsize>(header.size()));

    for (const auto &row : metrics)
      for (double value : row) {
        const auto bytes = std::bit_cast<std::array<char, sizeof(double)>>(value);
        if constexpr (std::endian::native == std::endian::little) {
          output.write(bytes.data(), bytes.size());
        } else {
          for (auto it = bytes.rbegin(); it != bytes.rend(); ++it)
            output.put(*it);
        }
      }
  }
};
