#include <chrono>
#include <iostream>
#include <ostream>
#include <sys/stat.h>
#include <sys/types.h>

#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>

#include "config.hpp"
#include "grid.hpp"
#include "particles.hpp"
#include "physical_converter.hpp"
#include "print_result.hpp"
#include "field_view_2d_result.hpp"

void run_3d() {
  std::cout << "Running 3D FDTD sim!!" << std::endl;

  const PhysicalConverter converter{1e-6};

  constexpr int nx = Config::nx;
  constexpr int ny = Config::ny;
  constexpr int nz = Config::nz;
  constexpr int print_interval = Config::print_interval;
  constexpr int slice_z = (nz - 1) / 2;
  constexpr double particle_x = (nx - 1) / 2.0;
  constexpr double particle_y = (ny - 1) / 8.0;
  constexpr double particle_z = (nz - 1) / 2.0;

  std::vector<SimulationResult<nx, ny, nz>> results;
  results.emplace_back(PrintResult{print_interval});
  results.emplace_back(
      FieldView2D<nx, ny, nz>("e-field", WhichField::E, Axis::Z, slice_z,
                              print_interval));

  const double sim_velocity = converter.to_sim_velocity(light_speed * 0.5);
  const double sim_charge = converter.to_sim_charge(-fundamental_charge);
  const double sim_mass = converter.to_sim_mass(electron_mass);

  std::cout << "v'=" << sim_velocity
    << ", q'=" << sim_charge
    << ", m'=" << sim_mass
    << std::endl;

  Sim<nx, ny, nz> grid(
    std::vector<Particle>{
      Particle{
        particle_x, particle_y, particle_z,
        0, sim_velocity, 0,
        sim_charge, sim_mass
      },
    },
    std::move(results)
  );

  grid.initialize();

  auto start = std::chrono::steady_clock::now();

  while (grid.time < Config::max_time)
    grid.step();
  
  grid.finish();

  auto finished = std::chrono::steady_clock::now();

  std::cout << "Sim finished!" << std::endl;
  std::cout << "took: " << std::chrono::duration_cast<std::chrono::milliseconds>(finished - start) << "ms" << std::endl;
}

int main() {
  run_3d();

  return 0;
}
