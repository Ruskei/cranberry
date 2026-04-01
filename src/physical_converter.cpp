#include "physical_converter.hpp"

/*
 * cell = grid_dx m => m = cell / grid_dx
 * timeunit = dt s => s = timeunit / dt
 * v m / s = v (cell / grid_dx) / (timeunit / dt) = v (cell / timeunit) (dt / grid_dx)
 */
double PhysicalConverter::to_sim_velocity(double physical_velocity) const {
  return physical_velocity * dt / grid_dx;
}
