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
#include "field_view_2d_result.hpp"
#include "grid.hpp"
#include "particle_distribution.hpp"
#include "particles.hpp"
#include "physical_converter.hpp"
#include "print_result.hpp"

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
  results.emplace_back(FieldView2D<nx, ny, nz>(
      "e-field", WhichField::E, Axis::Z, slice_z, print_interval));

  const double sim_velocity = converter.to_sim_velocity(light_speed * 0.5);
  const double sim_charge = converter.to_sim_charge(-fundamental_charge);
  const double sim_mass = converter.to_sim_mass(electron_mass);

  std::cout << "v'=" << sim_velocity << ", q'=" << sim_charge
            << ", m'=" << sim_mass << std::endl;

  std::vector<Particle> particles;

  double beam_lb_x = 14e-6, beam_lb_y = 8e-6, beam_lb_z = 14e-6;
  double beam_x = 4e-6, beam_y = 4e-6, beam_z = 4e-6;
  double beam_speed = light_speed * 0.9;

  double target_lb_x = 9e-6, target_lb_y = 24e-6, target_lb_z = 9e-6;
  double target_x = 14e-6, target_y = 20e-6, target_z = 14e-6;

  std::vector<Particle> beam{generate_particle_distribution(
      electron_species, 0, beam_speed, 0, 1, beam_lb_x, beam_lb_y,
      beam_lb_z, beam_lb_x + beam_x, beam_lb_y + beam_y, beam_lb_z + beam_z,
      1e22, converter)};

  std::vector<Particle> electrons{
    generate_particle_distribution(electron_species, 0, 0, 0, 1, target_lb_x, target_lb_y, target_lb_z,
        target_lb_x + target_x, target_lb_y + target_y, target_lb_z + target_z, 1e22, converter)};
  std::vector<Particle> positrons{
    generate_particle_distribution(positron_species, 0, 0, 0, 1, target_lb_x, target_lb_y, target_lb_z,
        target_lb_x + target_x, target_lb_y + target_y, target_lb_z + target_z, 1e22, converter)};

  particles.insert(particles.end(), electrons.begin(), electrons.end());
  particles.insert(particles.end(), positrons.begin(), positrons.end());
  particles.insert(particles.end(), beam.begin(), beam.end());

  Sim<nx, ny, nz> grid(particles, std::move(results));

  grid.initialize();

  auto start = std::chrono::steady_clock::now();

  grid.step();
  while (grid.time < Config::max_time)
    grid.step();

  grid.finish();

  auto finished = std::chrono::steady_clock::now();

  std::cout << "Sim finished!" << std::endl;
  std::cout << "took: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(finished -
                                                                     start)
            << "ms" << std::endl;
}

int main() {
  run_3d();

  return 0;
}
