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

#include "particle_view_3d_result.hpp"
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
  results.emplace_back(ParticleView3D{"particles", print_interval,
                                      ParticleVelocityOutput::Proper});

  std::vector<Particle> particles;

  const Vec3d target_lb{0, 0e-6, 0};
  const Vec3d target_ub{64e-6, 256e-6, 64e-6};
  const double target_density{1e23};

  const UniformDistribution target_electron_distribution{
      electron_species, target_lb, target_ub, target_density, 1};
  const UniformDistribution target_proton_distribution{
      proton_species, target_lb, target_ub, target_density, 1};

  const Vec3d beam_velocity{0, proper_velocity_to_coord_velocity(1'00 * light_speed), 0};
  std::cout << "beam_velocity=" << (beam_velocity.y / light_speed) << std::endl;
  const Vec3d beam_p{32e-6, 85e-6, 32e-6};
  const Vec3d beam_rms{2e-6, 4e-6, 2e-6};
  const double beam_q_total{-1e-10};
  const double beam_m_total{-beam_q_total / fundamental_charge * electron_mass};
  const int beam_num_particles{1000};

  const GaussianDistribution beam_distribution{
    beam_p, beam_rms,
    beam_q_total, beam_m_total,
    beam_num_particles
  };

  const Vec3d witness_velocity{beam_velocity};
  const Vec3d witness_p{32e-6, 20e-6, 32e-6};
  const Vec3d witness_rms{2e-6, 2e-6, 2e-6};
  const double witness_q_total{beam_q_total * 0.3};
  const double witness_m_total{-witness_q_total / fundamental_charge * electron_mass};
  const int witness_num_particles{1000};

  const GaussianDistribution witness_distribution{
    witness_p, witness_rms,
    witness_q_total, witness_m_total,
    witness_num_particles
  };

  std::vector<Particle> beam{generate_particle_distribution(beam_velocity, beam_distribution, converter)};
  std::vector<Particle> witness{generate_particle_distribution(witness_velocity, witness_distribution, converter)};

  std::vector<Particle> electrons{generate_particle_distribution(Vec3d{}, target_electron_distribution, converter)};
  std::vector<Particle> protons{generate_particle_distribution(Vec3d{}, target_proton_distribution, converter)};

  particles.insert(particles.end(), electrons.begin(), electrons.end());
  particles.insert(particles.end(), protons.begin(), protons.end());
  particles.insert(particles.end(), beam.begin(), beam.end());
  particles.insert(particles.end(), witness.begin(), witness.end());

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
