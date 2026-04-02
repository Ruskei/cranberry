#pragma once

#include <cmath>

inline const double light_speed{299'792'458};
inline const double electron_mass{9.109'383'7139e-31};
inline const double fundamental_charge{1.602e-19};
inline const double m0{electron_mass};
inline const double e0{8.8541878188e-12};

/*
 * okay let's say that a mass value of 1 in the simulation corresponds to the mass of an electron
 * a prime denotes a simulation value, so:
 * - m = M0 m'
 * - l = L0 l'
 * - t = T0 t'
 * - q = Q0 q'
 * - since c' = 1 = l'/t' => l' = t' => l / L0 = t / T0 => T0 = L0 t / l
 * 1 = [e0] = C^2 t^2 m^-1 l^-3
 * C^2 = e0 l^3 m t^-2 
 */

struct PhysicalConverter {
  double l0{1};
  double t0{1};
  double q0{1};

  PhysicalConverter(double l0) : l0{l0}, t0{l0 / light_speed}, q0{std::sqrt(e0 * m0 * l0 * light_speed * light_speed)} {}

  double to_sim_velocity(double physical_velocity) const;
  double to_sim_mass(double physical_mass) const;
  double to_sim_charge(double physical_charge) const;
};
