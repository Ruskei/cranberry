#include <cmath>
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

#include "abc.hpp"
#include "config.hpp"
#include "fdtd_types.hpp"
#include "grid.hpp"
#include "io.hpp"

struct Particle {
  double nx;
  double ny;
  double nz;

  double vx;
  double vy;
  double vz;

  double q;
  double m;
};

struct ParticleWriter {
  std::string name;

  ParticleWriter(std::string n) : name{n} {
    mkdir("out", 0755);
    mkdir(("out/" + name).c_str(), 0755);
  }

  void write_particles(const std::vector<Particle> &particles, double time) {
    vtkNew<vtkPoints> points;
    vtkNew<vtkFloatArray> velocity;
    velocity->SetName("velocity");
    velocity->SetNumberOfComponents(3);
    for (const auto &p : particles) {
      points->InsertNextPoint(p.nx, p.ny, p.nz);
      velocity->InsertNextTuple3(p.vx, p.vy, p.vz);
    }

    vtkNew<vtkPolyData> poly;
    poly->SetPoints(points);
    poly->GetPointData()->AddArray(velocity);
    poly->GetPointData()->SetActiveVectors("velocity");

    vtkNew<vtkXMLPolyDataWriter> writer;
    const std::string filename = "out/" + name + "/particles" + std::to_string(time) + ".vtp";
    writer->SetFileName(filename.c_str());
    writer->SetInputData(poly);
    writer->Write();
  }

  void write_all(double max_time, double dt) {
    std::ofstream output;
    output.open(("out/" + name + ".pvd").c_str());

    output << "<?xml version=\"1.0\"?>" << std::endl;
    output << "<VTKFile type=\"Collection\" version=\"0.1\" "
              "byte_order=\"LittleEndian\">"
           << std::endl;
    output << "  <Collection>" << std::endl;

    for (auto t{0.0}; t < max_time; t += dt) {
      output << "    <DataSet timestep=\"" << std::to_string(t)
             << "\" group=\"\" part=\"0\" file=\"" + name + "/particles"
             << std::to_string(t) << ".vtp\" name=\"particles\"/>" << std::endl;
    }

    output << "  </Collection>" << std::endl;
    output << "</VTKFile>" << std::endl;

    output.close();
  }
};

void test_point_io() {
  mkdir("out", 0755);

  constexpr int N = 100;

  vtkNew<vtkPoints> points;
  vtkNew<vtkFloatArray> velocity;

  velocity->SetName("velocity");
  velocity->SetNumberOfComponents(3);

  for (int i = 0; i < N; ++i) {
    double t = 0.1 * i;

    double x = std::cos(t);
    double y = std::sin(t);
    double z = 0.05 * i;

    points->InsertNextPoint(x, y, z);

    float vx = -std::sin(t);
    float vy = std::cos(t);
    float vz = 0.05f;

    velocity->InsertNextTuple3(vx, vy, vz);
  }

  vtkNew<vtkPolyData> poly;
  poly->SetPoints(points);
  poly->GetPointData()->AddArray(velocity);
  poly->GetPointData()->SetActiveVectors("velocity");

  vtkNew<vtkXMLPolyDataWriter> writer;
  writer->SetFileName("out/particles.vtp");
  writer->SetInputData(poly);
  writer->Write();
}

template <class Grid>
void push_particles(const Grid &grid, std::vector<Particle> &particles) {
  const double epsilon = 1e-11;
  const double dt = Config::dt;
  for (auto &p : particles) {
    const double nx = p.nx;
    const double ny = p.ny;
    const double nz = p.nz;

    double vx = p.vx;
    double vy = p.vy;
    double vz = p.vz;

    double vv = vx * vx + vy * vy + vz * vz;
    if (!std::isfinite(vv) || vv >= 1.0 - epsilon) {
      double scale = std::sqrt((1.0 - 1e-12) / std::max(vv, epsilon));
      vx *= scale;
      vy *= scale;
      vz *= scale;
      vv = vx * vx + vy * vy + vz * vz;
    }

    const double q = p.q;
    const double m = p.m;

    double lorentz = 1.0 / std::sqrt(1.0 - vv);

    const double ex = grid.ex.interpolate_at(nx, ny, nz);
    const double ey = grid.ey.interpolate_at(nx, ny, nz);
    const double ez = grid.ez.interpolate_at(nx, ny, nz);

    const double hx = grid.hx.interpolate_at(nx, ny, nz);
    const double hy = grid.hy.interpolate_at(nx, ny, nz);
    const double hz = grid.hz.interpolate_at(nx, ny, nz);

    const double ux = lorentz * vx;
    const double uy = lorentz * vy;
    const double uz = lorentz * vz;

    const double u_minus_x = ux + (q * dt / (2.0 * m)) * ex;
    const double u_minus_y = uy + (q * dt / (2.0 * m)) * ey;
    const double u_minus_z = uz + (q * dt / (2.0 * m)) * ez;
    const double umum =
        (u_minus_x * u_minus_x + u_minus_y * u_minus_y + u_minus_z * u_minus_z);

    const double lorentz_bar = 1.0 / sqrt(1.0 + umum);
    const double tx = q * dt / (2.0 * m) * hx * lorentz_bar;
    const double ty = q * dt / (2.0 * m) * hy * lorentz_bar;
    const double tz = q * dt / (2.0 * m) * hz * lorentz_bar;
    const double tt = tx * tx + ty * ty + tz * tz;

    const double u_prime_x = u_minus_x + (u_minus_y * tz - u_minus_z * ty);
    const double u_prime_y = u_minus_y + (u_minus_z * tx - u_minus_x * tz);
    const double u_prime_z = u_minus_z + (u_minus_x * ty - u_minus_y * tx);

    const double ax = 2 * tx / (1 + tt);
    const double ay = 2 * ty / (1 + tt);
    const double az = 2 * tz / (1 + tt);

    const double u_plus_x = u_minus_x + (u_prime_y * az - u_prime_z * ay);
    const double u_plus_y = u_minus_y + (u_prime_z * ax - u_prime_x * az);
    const double u_plus_z = u_minus_z + (u_prime_x * ay - u_prime_y * ax);

    const double u_next_x = u_plus_x + q * dt / (2 * m) * ex;
    const double u_next_y = u_plus_y + q * dt / (2 * m) * ey;
    const double u_next_z = u_plus_z + q * dt / (2 * m) * ez;
    const double unun =
        (u_next_x * u_next_x + u_next_y * u_next_y + u_next_z * u_next_z);

    const double lorentz_next = std::sqrt(1.0 + unun);
    p.vx = u_next_x / lorentz_next;
    p.vy = u_next_y / lorentz_next;
    p.vz = u_next_z / lorentz_next;

    p.nx += p.vx * dt;
    p.ny += p.vy * dt;
    p.nz += p.vz * dt;
  }
}

void run_3d() {
  test_point_io();

  std::vector<Particle> particles(1);
  particles[0] = Particle{10, 20, 20, 0, 0, 0, 1, 1};

  ParticleWriter particle_writer("particle_pusher");
  particle_writer.write_particles(particles, -1);

  std::cout << "Running 3D FDTD sim!!" << std::endl;

  const int size = Config::size;

  // const Writer2D<size> writer("2d_animation", (size - 1) / 2);
  const Writer3D<size> writer;

  ABC<size> abc;
  Grid<size> grid;

  grid.setup_coefficients();

  const double dt = Config::dt;

  for (auto time{0.0}; time < Config::max_time; time += dt) {
    grid.update_h();
    grid.update_e();
    abc.apply(grid);

    push_particles(grid, particles);
    particle_writer.write_particles(particles, time);

    grid.update_h();
    grid.update_e();
    abc.apply(grid);

    double tau = time - Config::t0;
    double a = M_PI * M_PI * Config::f0 * Config::f0 * tau * tau;
    grid.ex(grid.ex.nx() / 2, grid.ex.ny() / 2, grid.ex.nz() / 2) +=
        Config::A * (1.0 - 2.0 * a) * std::exp(-a);

    writer.write_timestep(grid, time);
  }

  writer.write_all(Config::max_time, Config::dt);
  particle_writer.write_all(Config::max_time, Config::dt);

  std::cout << "Finished 3D simulation!" << std::endl;
}

int main() {
  run_3d();

  return 0;
}
