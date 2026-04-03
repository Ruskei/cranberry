#include <cassert>
#include <vector>

#include "particle_distribution.hpp"
#include "particles.hpp"
#include "physical_converter.hpp"

std::vector<Particle> generate_particle_distribution(
    ParticleSpecies species, double vx, double vy, double vz,
    double macroparticles_per_cell,
    double lower_bound_x, double lower_bound_y, double lower_bound_z,
    double upper_bound_x, double upper_bound_y, double upper_bound_z,
    double density,
    PhysicalConverter converter
    ) {
  double sim_vx = converter.to_sim_velocity(vx);
  double sim_vy = converter.to_sim_velocity(vy);
  double sim_vz = converter.to_sim_velocity(vz);

  double range_x = upper_bound_x - lower_bound_x;
  double range_y = upper_bound_y - lower_bound_y;
  double range_z = upper_bound_z - lower_bound_z;
  assert(range_x > 0 && range_y > 0 && range_z > 0);

  double volume{range_x * range_y * range_z};
  double particles{density * volume};
  int macroparticles{static_cast<int>(volume / converter.l0 / converter.l0 / converter.l0)};
  double physical_particles_per_macroparticle = particles / macroparticles;

  std::vector<Particle> result(macroparticles);

  double spacing = converter.l0 / macroparticles_per_cell;
  for (double rx{lower_bound_x}; rx < upper_bound_x; rx += spacing)
    for (double ry{lower_bound_y}; ry < upper_bound_y; ry += spacing)
      for (double rz{lower_bound_z}; rz < upper_bound_z; rz += spacing) {
        result.emplace_back(Particle{
          rx / converter.l0, ry / converter.l0, rz / converter.l0,
          sim_vx, sim_vy, sim_vz,
          converter.to_sim_charge(species.q) * physical_particles_per_macroparticle,
          converter.to_sim_mass(species.m) * physical_particles_per_macroparticle
        });
      }

  return result;
}
