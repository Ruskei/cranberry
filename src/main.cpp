#include <chrono>
#include <iostream>
#include <ostream>
#include <ratio>
#include <sys/stat.h>
#include <sys/types.h>

#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>

#include "abc.hpp"
#include "config.hpp"
#include "grid.hpp"
#include "io.hpp"
#include "particles.hpp"
#include "progress_bar.hpp"

void run_3d() {
  std::cout << "Running 3D FDTD sim!!" << std::endl;

  const int size = Config::size;

  ABC<size> abc;
  Grid<size> grid(std::vector<Particle>{
      Particle{64, 50, 64, 0, 0.5, 0, 1, 1},
      Particle{64, 70, 64, 0, -0.5, 0, -1, 1},
  });

  const Writer2D<size> writer("2d_animation", (size - 1) / 2);
  // const Writer3D<size> writer;
  ParticleWriter particle_writer("particle_pusher");
  particle_writer.write_particles(grid.particles, -1);

  grid.setup_coefficients();

  const double dt = Config::dt;

  /*
   * half-update E^(n    ) -> E^(n+1/2) with H^(n    )
   * half-update H^(n    ) -> H^(n+1/2) with E^(n+1/2)
   * PIC
   * half-update H^(n+1/2) -> H^(n+1  ) with E^(n+1/2)
   * half-update E^(n+1/2) -> E^(n    ) with H^(n+1/2)
   */
  for (auto time{0.0}; time < Config::max_time; time += dt) {
    grid.half_update_e();
    grid.half_update_h();
    abc.apply(grid);

    grid.step_particles();
    grid.push_particles();
    grid.deposit_currents();
    grid.check_currents();
    grid.smooth_currents();

    // particle_writer.write_particles(grid.particles, time);

    grid.half_update_h();
    grid.half_update_e();
    abc.apply(grid);

    grid.check_gauss();

    // if (time > 30.0) {
    // double tau = time - Config::t0;
    // double a = M_PI * M_PI * Config::f0 * Config::f0 * tau * tau;
    // grid.ez(grid.ez.nx() / 2, grid.ez.ny() / 2, grid.ez.nz() / 2) +=
    // Config::A * (1.0 - 2.0 * a) * std::exp(-a);
    // }

    writer.write_timestep(grid, time);

    if (static_cast<int>(time / dt) % Config::print_interval == 0)
      print_progress(time / Config::max_time, 20);
  }

  writer.write_all(Config::max_time, Config::dt);
  particle_writer.write_all(Config::max_time, Config::dt);

  print_sim_finished();
}

int main() {
  run_3d();

  return 0;
}
