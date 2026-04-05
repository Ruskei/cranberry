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

#include "linalg.hpp"
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
  constexpr int slice_z = (nx - 1) / 2;

  std::vector<SimulationResult<nx, ny, nz>> results;
  results.emplace_back(PrintResult{print_interval});
  results.emplace_back(FieldView2D<nx, ny, nz>(
      "e-field", WhichField::E, Axis::Z, slice_z, print_interval));

  std::vector<Particle> particles;

  const Vec3d target_lb{8e-6, 16e-6, 8e-6};
  const Vec3d target_ub{56e-6, 124e-6, 56e-6};
  const double target_density{1e23};

  const UniformDistribution target_electron_distribution{
      electron_species, target_lb, target_ub, target_density, 1};
  const UniformDistribution target_proton_distribution{
      proton_species, target_lb, target_ub, target_density, 1};

  const Vec3d beam_velocity{0, light_speed * 0.9, 0};
  const Vec3d beam_p{32e-6, 12e-6, 32e-6};
  const Vec3d beam_rms{2e-6, 4e-6, 2e-6};
  const double beam_q_total{-1e-10};
  const double beam_m_total{1e10};
  const int beam_num_particles{1000};

  const GaussianDistribution beam_distribution{
    beam_p, beam_rms,
    beam_q_total, beam_m_total,
    beam_num_particles
  };

  std::vector<Particle> beam{generate_particle_distribution(beam_velocity, beam_distribution, converter)};

  std::vector<Particle> electrons{generate_particle_distribution(Vec3d{}, target_electron_distribution, converter)};
  std::vector<Particle> positrons{generate_particle_distribution(Vec3d{}, target_proton_distribution, converter)};

  particles.insert(particles.end(), electrons.begin(), electrons.end());
  particles.insert(particles.end(), positrons.begin(), positrons.end());
  particles.insert(particles.end(), beam.begin(), beam.end());

  Sim<nx, ny, nz> grid(particles, std::move(results));

  grid.initialize();

  auto start = std::chrono::steady_clock::now();

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
