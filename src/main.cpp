#include <cmath>
#include <iostream>
#include <ostream>
#include <sys/stat.h>
#include <sys/types.h>

#include "config.hpp"
#include "fdtd_types.hpp"
#include "grid.hpp"
#include "io.hpp"
#include "abc.hpp"

void run_3d() {
  std::cout << "Running 3D FDTD sim!!" << std::endl;

  const int size = Config::size;

  const Writer2D<size> writer((size - 1) / 2);

  ABC<size> abc;
  Grid<size> grid;

  grid.setup_coefficients();

  for (auto time{0}; time < Config::max_time; ++time) {
    const double t = time * Config::courants;
    grid.update_h();
    grid.update_e();

    double tau = t - Config::t0;
    double a = M_PI * M_PI * Config::f0 * Config::f0 * tau * tau;
    grid.ex(grid.ex.nx() / 2, grid.ex.ny() / 2, grid.ex.nz() / 2) +=
        (1.0 - 2.0 * a) * std::exp(-a);

    abc.apply(grid);

    writer.write_timestep(grid, time);
  }

  writer.write_all(Config::max_time);

  std::cout << "Finished 3D simulation!" << std::endl;
}

int main() {
  run_3d();

  return 0;
}
