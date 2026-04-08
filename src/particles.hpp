#pragma once

#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>
#include <cmath>

namespace {
  inline constexpr int particle_radius{2};
}

struct Particle {
  double px;
  double py;
  double pz;

  double p_prev_x;
  double p_prev_y;
  double p_prev_z;

  double vx;
  double vy;
  double vz;

  double ux;
  double uy;
  double uz;

  double q;
  double m;

  int age{};

  Particle() = default;

  Particle(double px, double py, double pz, double vx, double vy, double vz,
           double q, double m)
      : px{px}, py{py}, pz{pz}, p_prev_x{px}, p_prev_y{py}, p_prev_z{pz},
        vx{vx}, vy{vy}, vz{vz}, ux{}, uy{}, uz{}, q{q}, m{m} {
    const double gamma =
        1.0 / std::sqrt(1.0 - (vx * vx + vy * vy + vz * vz));
    ux = gamma * vx;
    uy = gamma * vy;
    uz = gamma * vz;
  }
};

enum class CartesianComponent { x, y, z };

struct ParticleWriter {
  std::string name;

  ParticleWriter(std::string n) : name{n} {
    mkdir("out", 0755);
    mkdir(("out/" + name).c_str(), 0755);
  }

  void write_particles(const std::vector<Particle> &particles, double time);

  void write_all(double max_time, double dt);
};

double form_factor(double px, double nx);

double form_factor_diff_helper(double s_o_x, double s_o_y, double s_o_z,
                               double s_n_x, double s_n_y, double s_n_z);

// template <CartesianComponent C>
// double form_factor_diff(const Particle &particle,
//     double a, double b,
//     double c, double d,
//     double e, double f
//   ) {
//   const double px = particle.px;
//   const double py = particle.py;
//   const double pz = particle.pz;
//
//   const double ppx = particle.p_prev_x;
//   const double ppy = particle.p_prev_y;
//   const double ppz = particle.p_prev_z;
//
//   if constexpr (C == CartesianComponent::x) {
//     return form_factor_diff_helper(form_factor(ppx, i), form_factor(ppy, j),
//                                    form_factor(ppz, k), form_factor(px, i),
//                                    form_factor(py, j), form_factor(pz, k));
//   } else if constexpr (C == CartesianComponent::y) {
//     return form_factor_diff_helper(form_factor(ppy, j), form_factor(ppx, i),
//                                    form_factor(ppz, k), form_factor(py, j),
//                                    form_factor(px, i), form_factor(pz, k));
//   } else {
//     return form_factor_diff_helper(form_factor(ppz, k), form_factor(ppy, j),
//                                    form_factor(ppx, i), form_factor(pz, k),
//                                    form_factor(py, j), form_factor(px, i));
//   }
// }

double form_factor_x1(const Particle &particle, double idx);
double form_factor_x2(const Particle &particle, double idx);
double form_factor_y1(const Particle &particle, double idx);
double form_factor_y2(const Particle &particle, double idx);
double form_factor_z1(const Particle &particle, double idx);
double form_factor_z2(const Particle &particle, double idx);

template <CartesianComponent C>
double form_factor_diff(
    double x1, double x2,
    double y1, double y2,
    double z1, double z2
  ) {
  if constexpr (C == CartesianComponent::x) {
    return form_factor_diff_helper(x2, y2,
                                   z2, x1,
                                   y1, z1);
  } else if constexpr (C == CartesianComponent::y) {
    return form_factor_diff_helper(y2, x2,
                                   z2, y1,
                                   x1, z1);
  } else {
    return form_factor_diff_helper(z2, y2,
                                   x2, z1,
                                   y1, x1);
  }
}

template <CartesianComponent C>
double form_factor_diff(const Particle &particle, double i, double j, double k) {
  const double px = particle.px;
  const double py = particle.py;
  const double pz = particle.pz;

  const double ppx = particle.p_prev_x;
  const double ppy = particle.p_prev_y;
  const double ppz = particle.p_prev_z;

  if constexpr (C == CartesianComponent::x) {
    return form_factor_diff_helper(form_factor(ppx, i), form_factor(ppy, j),
                                   form_factor(ppz, k), form_factor(px, i),
                                   form_factor(py, j), form_factor(pz, k));
  } else if constexpr (C == CartesianComponent::y) {
    return form_factor_diff_helper(form_factor(ppy, j), form_factor(ppx, i),
                                   form_factor(ppz, k), form_factor(py, j),
                                   form_factor(px, i), form_factor(pz, k));
  } else {
    return form_factor_diff_helper(form_factor(ppz, k), form_factor(ppy, j),
                                   form_factor(ppx, i), form_factor(pz, k),
                                   form_factor(py, j), form_factor(px, i));
  }
}
