#pragma once

#include "linalg.hpp"
#include "particles.hpp"
#include "physical_converter.hpp"

struct ParticleSpecies {
  double q, m;
};

struct UniformDistribution {
  ParticleSpecies species;
  Vec3d lower_bound;
  Vec3d upper_bound;
  double density;
  double macroparticles_per_cell;
};

struct GaussianDistribution {
  Vec3d mean_pos;
  Vec3d rms;
  double q_total;
  double m_total;
  int num_particles;
};

using Distribution = std::variant<UniformDistribution, GaussianDistribution>;

inline const ParticleSpecies electron_species{-fundamental_charge, electron_mass};
inline const ParticleSpecies positron_species{fundamental_charge, electron_mass};

inline const ParticleSpecies proton_species{fundamental_charge, proton_mass};

std::vector<Particle> generate_particle_distribution(
    Vec3d velocity,
    Distribution distribution,
    PhysicalConverter converter
);
