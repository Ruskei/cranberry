#pragma once

#include "particles.hpp"
#include "physical_converter.hpp"

struct ParticleSpecies {
  double q, m;
};

inline const ParticleSpecies electron_species{-fundamental_charge, electron_mass};

std::vector<Particle> generate_particle_distribution(
    ParticleSpecies species, double macroparticles_per_cell,
    double lower_bound_x, double lower_bound_y, double lower_bound_z,
    double upper_bound_x, double upper_bound_y, double upper_bound_z,
    double density,
    PhysicalConverter converter
);
