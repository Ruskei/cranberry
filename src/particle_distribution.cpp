#include <cassert>
#include <vector>
#include <random>

#include "particle_distribution.hpp"
#include "particles.hpp"
#include "physical_converter.hpp"

template<class... Ts>
struct overloads : Ts... { using Ts::operator()...; };

std::vector<Particle> generate_particle_distribution(
    Vec3d velocity,
    Distribution distribution,
    PhysicalConverter converter
) {
  Vec3d sim_v = Vec3d{
    converter.to_sim_velocity(velocity.x),
    converter.to_sim_velocity(velocity.y),
    converter.to_sim_velocity(velocity.z),
  };

  const auto visitor = overloads(
    [&](UniformDistribution &distribution) -> std::vector<Particle> {
      const double density{distribution.density};
      const Vec3d lb{distribution.lower_bound};
      const Vec3d ub{distribution.upper_bound};
      const double n{distribution.macroparticles_per_cell};
      const ParticleSpecies species{distribution.species};

      Vec3d range = distribution.upper_bound - distribution.lower_bound;
      
      double volume{range.x * range.y * range.z};
      double particles{density * volume};
      int macroparticles{static_cast<int>(volume / converter.l0 / converter.l0 / converter.l0)};
      double physical_particles_per_macroparticle = particles / macroparticles;

      std::vector<Particle> result;

      double spacing = converter.l0 / n;
      for (double rx{lb.x}; rx < ub.x; rx += spacing)
        for (double ry{lb.y}; ry < ub.y; ry += spacing)
          for (double rz{lb.z}; rz < ub.z; rz += spacing) {
            result.emplace_back(Particle{
              rx / converter.l0, ry / converter.l0, rz / converter.l0,
              sim_v.x, sim_v.y, sim_v.z,
              converter.to_sim_charge(species.q) * physical_particles_per_macroparticle,
              converter.to_sim_mass(species.m) * physical_particles_per_macroparticle
            });
          }

      return result;
    },
    [&](GaussianDistribution &distribution) -> std::vector<Particle> {
      const Vec3d p = distribution.mean_pos;
      const Vec3d rms = distribution.rms;
      const int n = distribution.num_particles;

      const double q_t{distribution.q_total};
      const double q_per{q_t / n};
      const double p_q_per{converter.to_sim_charge(q_per)};

      const double m_t{distribution.m_total};
      const double m_per{m_t / n};
      const double p_m_per{converter.to_sim_mass(m_per)};

      auto x_dist = std::normal_distribution{p.x, rms.x};
      auto y_dist = std::normal_distribution{p.y, rms.y};
      auto z_dist = std::normal_distribution{p.z, rms.z};

      const double l0 = converter.l0;
      
      std::random_device rd{};
      std::mt19937 gen{rd()};

      std::vector<Particle> result;

      for (int i{0}; i < n; ++i)
        result.emplace_back(Particle{
          x_dist(gen) / l0, y_dist(gen) / l0, z_dist(gen) / l0,
          sim_v.x, sim_v.y, sim_v.z,
          p_q_per, p_m_per
        });
      
      return result;
    }
  );

  return std::visit(visitor, distribution);
}
