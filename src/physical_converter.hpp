#pragma once

inline const double light_speed = 299'792'458;

struct PhysicalConverter {
  double grid_dx{1};
  double dt{1};

  PhysicalConverter(double grid_dx) : grid_dx{grid_dx}, dt{grid_dx / light_speed} {};

  double to_sim_velocity(double physical_velocity) const;
};
