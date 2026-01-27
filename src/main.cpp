#include <algorithm>
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
#include "grid.hpp"
#include "io.hpp"
#include "progress_bar.hpp"

struct Particle {
  double px;
  double py;
  double pz;

  double p_prev_x;
  double p_prev_y;
  double p_prev_z;

  double vx;
  double vy;
  double vz;

  double q;
  double m;

  double weight{1.0};
  double weight_prev{1.0};

  int age{};

  Particle() = default;

  Particle(double px, double py, double pz, double vx, double vy, double vz,
           double q, double m)
      : px{px}, py{py}, pz{pz}, p_prev_x{px}, p_prev_y{py}, p_prev_z{pz},
        vx{vx}, vy{vy}, vz{vz}, q{q}, m{m} {}
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
      points->InsertNextPoint(p.px, p.py, p.pz);
      velocity->InsertNextTuple3(p.vx, p.vy, p.vz);
    }

    vtkNew<vtkPolyData> poly;
    poly->SetPoints(points);
    poly->GetPointData()->AddArray(velocity);
    poly->GetPointData()->SetActiveVectors("velocity");

    vtkNew<vtkXMLPolyDataWriter> writer;
    const std::string filename =
        "out/" + name + "/particles" + std::to_string(time) + ".vtp";
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

void push_particles(std::vector<Particle> &particles) {
  const double dt = Config::dt;
  for (auto &p : particles) {
    // const double nx = p.px;
    // const double ny = p.py;
    // const double nz = p.pz;
    //
    // double vx = p.vx;
    // double vy = p.vy;
    // double vz = p.vz;
    //
    // double vv = vx * vx + vy * vy + vz * vz;
    // if (!std::isfinite(vv) || vv >= 1.0 - epsilon) {
    //   double scale = std::sqrt((1.0 - 1e-12) / std::max(vv, epsilon));
    //   vx *= scale;
    //   vy *= scale;
    //   vz *= scale;
    //   vv = vx * vx + vy * vy + vz * vz;
    // }
    //
    // const double q = p.q;
    // const double m = p.m;
    //
    // double lorentz = 1.0 / std::sqrt(1.0 - vv);
    //
    // const double ex = grid.ex.interpolate_at(nx, ny, nz);
    // const double ey = grid.ey.interpolate_at(nx, ny, nz);
    // const double ez = grid.ez.interpolate_at(nx, ny, nz);
    //
    // const double hx = grid.hx.interpolate_at(nx, ny, nz);
    // const double hy = grid.hy.interpolate_at(nx, ny, nz);
    // const double hz = grid.hz.interpolate_at(nx, ny, nz);
    //
    // const double ux = lorentz * vx;
    // const double uy = lorentz * vy;
    // const double uz = lorentz * vz;
    //
    // const double u_minus_x = ux + (q * dt / (2.0 * m)) * ex;
    // const double u_minus_y = uy + (q * dt / (2.0 * m)) * ey;
    // const double u_minus_z = uz + (q * dt / (2.0 * m)) * ez;
    // const double umum =
    //     (u_minus_x * u_minus_x + u_minus_y * u_minus_y + u_minus_z *
    //     u_minus_z);
    //
    // const double lorentz_bar = 1.0 / sqrt(1.0 + umum);
    // const double tx = q * dt / (2.0 * m) * hx * lorentz_bar;
    // const double ty = q * dt / (2.0 * m) * hy * lorentz_bar;
    // const double tz = q * dt / (2.0 * m) * hz * lorentz_bar;
    // const double tt = tx * tx + ty * ty + tz * tz;
    //
    // const double u_prime_x = u_minus_x + (u_minus_y * tz - u_minus_z * ty);
    // const double u_prime_y = u_minus_y + (u_minus_z * tx - u_minus_x * tz);
    // const double u_prime_z = u_minus_z + (u_minus_x * ty - u_minus_y * tx);
    //
    // const double ax = 2 * tx / (1 + tt);
    // const double ay = 2 * ty / (1 + tt);
    // const double az = 2 * tz / (1 + tt);
    //
    // const double u_plus_x = u_minus_x + (u_prime_y * az - u_prime_z * ay);
    // const double u_plus_y = u_minus_y + (u_prime_z * ax - u_prime_x * az);
    // const double u_plus_z = u_minus_z + (u_prime_x * ay - u_prime_y * ax);
    //
    // const double u_next_x = u_plus_x + q * dt / (2 * m) * ex;
    // const double u_next_y = u_plus_y + q * dt / (2 * m) * ey;
    // const double u_next_z = u_plus_z + q * dt / (2 * m) * ez;
    // const double unun =
    //     (u_next_x * u_next_x + u_next_y * u_next_y + u_next_z * u_next_z);
    //
    // const double lorentz_next = std::sqrt(1.0 + unun);
    // p.vx = u_next_x / lorentz_next;
    // p.vy = u_next_y / lorentz_next;
    // p.vz = u_next_z / lorentz_next;
    //
    // p.p_prev_x = p.px;
    // p.p_prev_y = p.py;
    // p.p_prev_z = p.pz;
    //
    // p.px += p.vx * dt;
    // p.py += p.vy * dt;
    // p.pz += p.vz * dt;

    p.p_prev_x = p.px;
    p.p_prev_y = p.py;
    p.p_prev_z = p.pz;

    p.px += p.vx * dt;
    p.py += p.vy * dt;
    p.pz += p.vz * dt;
  }
}

double form_factor(double px, int nx) {
  const int i = static_cast<int>(std::floor(px));

  if (nx < i || nx > i + 1)
    return 0.0;

  const double fx = px - i;

  return (i == nx) ? 1.0 - fx : fx;
}

double form_factor_diff_helper(double s_o_x, double s_o_y, double s_o_z,
                               double s_n_x, double s_n_y, double s_n_z) {
  return (s_n_x * s_n_y * s_n_z - s_o_x * s_n_y * s_n_z +
          s_n_x * s_o_y * s_o_z - s_o_x * s_o_y * s_o_z) /
             3.0 +
         (s_n_x * s_o_y * s_n_z - s_o_x * s_o_y * s_n_z +
          s_n_x * s_n_y * s_o_z - s_o_x * s_n_y * s_o_z) /
             6.0;
}

enum class CartesianComponent { x, y, z };

template <CartesianComponent C>
double form_factor_diff(const Particle &particle, int i, int j, int k) {
  const double px = particle.px;
  const double py = particle.py;
  const double pz = particle.pz;

  const double ppx = particle.p_prev_x;
  const double ppy = particle.p_prev_y;
  const double ppz = particle.p_prev_z;

  if constexpr (C == CartesianComponent::x) {
    return form_factor_diff_helper(form_factor(ppx, i), form_factor(ppy, j),
                                   form_factor(ppz, k), form_factor(px, i),
                                   form_factor(py, j), form_factor(pz, k));
  } else if constexpr (C == CartesianComponent::y) {
    return form_factor_diff_helper(form_factor(ppy, j), form_factor(ppx, i),
                                   form_factor(ppz, k), form_factor(py, j),
                                   form_factor(px, i), form_factor(pz, k));
  } else {
    return form_factor_diff_helper(form_factor(ppz, k), form_factor(ppy, j),
                                   form_factor(ppx, i), form_factor(pz, k),
                                   form_factor(py, j), form_factor(px, i));
  }
}

template <class Grid>
void deposit_currents(Grid &grid, std::vector<Particle> &particles) {
  std::fill(grid.jx.v.begin(), grid.jx.v.end(), 0);
  std::fill(grid.jy.v.begin(), grid.jy.v.end(), 0);
  std::fill(grid.jz.v.begin(), grid.jz.v.end(), 0);

  for (const auto &p : particles) {
    const double px = p.px;
    const double py = p.py;
    const double pz = p.pz;

    const double ppx = p.p_prev_x;
    const double ppy = p.p_prev_y;
    const double ppz = p.p_prev_z;

    const int i = static_cast<int>(std::floor(px));
    const int j = static_cast<int>(std::floor(py));
    const int k = static_cast<int>(std::floor(pz));

    const int pi = static_cast<int>(std::floor(ppx));
    const int pj = static_cast<int>(std::floor(ppy));
    const int pk = static_cast<int>(std::floor(ppz));

    const int is = std::max(0, std::min(i, pi));
    const int js = std::max(0, std::min(j, pj));
    const int ks = std::max(0, std::min(k, pk));

    const int ie = std::min(grid.jx.nx() - 1, std::max(i, pi) + 1);
    const int je = std::min(grid.jy.ny() - 1, std::max(j, pj) + 1);
    const int ke = std::min(grid.jz.nz() - 1, std::max(k, pk) + 1);

    const double c = -p.q / Config::courants;
    const double move_co = 0.5 * (p.weight + p.weight_prev) * c;
    const double src_co = (p.weight - p.weight_prev) * -p.q;

    for (auto y{js}; y <= je; ++y)
      for (auto z{ks}; z <= ke; ++z) {
        double move_sum = 0.0;
        double src_sum = 0.0;
        for (auto x{is}; x <= ie; ++x) {
          move_sum +=
              move_co * form_factor_diff<CartesianComponent::x>(p, x, y, z);
          const double sx =
              0.5 * (form_factor(p.p_prev_x, x) + form_factor(p.px, x));
          const double sy =
              0.5 * (form_factor(p.p_prev_y, y) + form_factor(p.py, y));
          const double sz =
              0.5 * (form_factor(p.p_prev_z, z) + form_factor(p.pz, z));
          src_sum += src_co * sx * sy * sz / 3.0;

          grid.jx(x, y, z) += move_sum + src_sum;
        }
      }

    for (auto x{is}; x <= ie; ++x)
      for (auto z{ks}; z <= ke; ++z) {
        double move_sum = 0.0;
        double src_sum = 0.0;
        for (auto y{js}; y <= je; ++y) {
          move_sum +=
              move_co * form_factor_diff<CartesianComponent::y>(p, x, y, z);
          const double sx =
              0.5 * (form_factor(p.p_prev_x, x) + form_factor(p.px, x));
          const double sy =
              0.5 * (form_factor(p.p_prev_y, y) + form_factor(p.py, y));
          const double sz =
              0.5 * (form_factor(p.p_prev_z, z) + form_factor(p.pz, z));
          src_sum += src_co * sx * sy * sz / 3.0;

          grid.jy(x, y, z) += move_sum + src_sum;
        }
      }

    for (auto x{is}; x <= ie; ++x)
      for (auto y{js}; y <= je; ++y) {
        double move_sum = 0.0;
        double src_sum = 0.0;
        for (auto z{ks}; z <= ke; ++z) {
          move_sum +=
              move_co * form_factor_diff<CartesianComponent::z>(p, x, y, z);
          const double sx =
              0.5 * (form_factor(p.p_prev_x, x) + form_factor(p.px, x));
          const double sy =
              0.5 * (form_factor(p.p_prev_y, y) + form_factor(p.py, y));
          const double sz =
              0.5 * (form_factor(p.p_prev_z, z) + form_factor(p.pz, z));
          src_sum += src_co * sx * sy * sz / 3.0;

          grid.jz(x, y, z) += move_sum + src_sum;
        }
      }
  }
}

template <class Grid>
void check_currents(Grid &grid, std::vector<Particle> &particles) {
  for (const auto &p : particles) {
    double residual{};

    const double px = p.px;
    const double py = p.py;
    const double pz = p.pz;

    const double ppx = p.p_prev_x;
    const double ppy = p.p_prev_y;
    const double ppz = p.p_prev_z;

    const int i = static_cast<int>(std::floor(px));
    const int j = static_cast<int>(std::floor(py));
    const int k = static_cast<int>(std::floor(pz));

    const int pi = static_cast<int>(std::floor(ppx));
    const int pj = static_cast<int>(std::floor(ppy));
    const int pk = static_cast<int>(std::floor(ppz));

    const int is = std::max(0, std::min(i, pi));
    const int js = std::max(0, std::min(j, pj));
    const int ks = std::max(0, std::min(k, pk));

    const int ie = std::min(grid.jx.nx() - 1, std::max(i, pi) + 1);
    const int je = std::min(grid.jy.ny() - 1, std::max(j, pj) + 1);
    const int ke = std::min(grid.jz.nz() - 1, std::max(k, pk) + 1);

    for (auto x{is}; x <= ie; ++x)
      for (auto y{js}; y <= je; ++y)
        for (auto z{ks}; z <= ke; ++z) {
          const double charge = p.q * form_factor(px, x) * form_factor(py, y) *
                                form_factor(pz, z);
          const double charge_prev = p.q * form_factor(ppx, x) *
                                     form_factor(ppy, y) * form_factor(ppz, z);

          const double clamped_x = std::max(0, x - 1);
          const double clamped_y = std::max(0, y - 1);
          const double clamped_z = std::max(0, z - 1);

          const double div_J = (grid.jx(x, y, z) - grid.jx(clamped_x, y, z)) +
                               (grid.jy(x, y, z) - grid.jy(x, clamped_y, z)) +
                               (grid.jz(x, y, z) - grid.jz(x, y, clamped_z));

          // if (std::abs(div_J) > 0.0)
          //   std::cout << "∇∙J_(" << x << "," << y << "," << z << ")=" <<
          //   div_J << std::endl;

          residual += ((charge - charge_prev) / Config::dt + div_J);
        }

    if (std::abs(residual) > 1e-11)
      std::cout << "continuity residual: " << residual << " at age: " << p.age
                << std::endl;
  }
}

void step_particles(std::vector<Particle> &particles, double time) {
  if (std::abs(time / Config::dt - 3.0) < 1e-11) {
    Particle p = Particle{50, 10, 50, 0, 0.55, 0, 1, 1};
    // p.weight = 1.0;
    // p.weight_prev = 0.0;
    particles.push_back(p);
  }

  for (auto it = particles.begin(); it != particles.end();) {
    // if (it->age > 50.0) {
    //   it = particles.erase(it);
    //   continue;
    // }

    auto &p = *it;

    switch (p.age) {
    // case 1:
    // p.weight = 1.0;
    // p.weight_prev = 1.0;
    // break;
    case 5:
      p.vy = 0.5;
      break;
    case 45:
      p.vy = 0.0;
      break;
      // case 50:
      // p.weight = 0.0;
      // p.weight_prev = 1.0;
      // break;
    }

    p.age++;
    ++it;
  }
}

void run_3d() {
  test_point_io();

  std::vector<Particle> particles(0);

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
    grid.half_update_e();
    grid.half_update_h();
    abc.apply(grid);

    step_particles(particles, time);
    push_particles(particles);
    deposit_currents(grid, particles);
    check_currents(grid, particles);

    particle_writer.write_particles(particles, time);

    grid.half_update_e();
    grid.half_update_h();
    abc.apply(grid);

    // if (time > 30.0) {
    // double tau = time - Config::t0;
    // double a = M_PI * M_PI * Config::f0 * Config::f0 * tau * tau;
    // grid.ez(grid.ez.nx() / 2, grid.ez.ny() / 2, grid.ez.nz() / 2) +=
    //     Config::A * (1.0 - 2.0 * a) * std::exp(-a);
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
