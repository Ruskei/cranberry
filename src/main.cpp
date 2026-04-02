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

#include "particle_distribution.hpp"
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
  constexpr int slice_z = 15;

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

  std::vector<Particle> particles;

  std::vector<Particle> electrons{generate_particle_distribution(
      electron_species, 0.2, 
      10e-6, 32e-6, 10e-6,
      20e-6, 64e-6, 20e-6, 
      1e22, converter
  )};

  particles.insert(particles.end(), electrons.begin(), electrons.end());

  Sim<nx, ny, nz> grid(
    particles,
    std::move(results)
  );

  grid.initialize();

  auto start = std::chrono::steady_clock::now();

  grid.step();
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
