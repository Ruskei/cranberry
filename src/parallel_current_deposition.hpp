#pragma once

#include <condition_variable>
#include <mutex>
#include <thread>

#include "config.hpp"
#include "fdtd_types.hpp"
#include "particles.hpp"

template <int NX, int NY, int NZ> struct ParallelCurrentDepositor {
  struct Worker {
    std::thread thread;
    int start{}, end{};
    JField<NX, NY, NZ> J;
  };

  void worker_loop(int thread_id) {
    int local_job_id{};

    for (;;) {
      {
        std::unique_lock<std::mutex> lock(m);
        state_cv.wait(lock, [&] { return stop || local_job_id != job_id; });
        if (stop)
          return;

        local_job_id = job_id;
      }

      if (thread_id < active_workers) {
        JField<NX, NY, NZ> &J = workers[thread_id].J;

        std::fill(J.x.v.begin(), J.x.v.end(), 0);
        std::fill(J.y.v.begin(), J.y.v.end(), 0);
        std::fill(J.z.v.begin(), J.z.v.end(), 0);

        constexpr int max_points = particle_radius * 2 + 1;
        std::array<double, max_points> x1_form_factors{};
        std::array<double, max_points> x2_form_factors{};
        std::array<double, max_points> y1_form_factors{};
        std::array<double, max_points> y2_form_factors{};
        std::array<double, max_points> z1_form_factors{};
        std::array<double, max_points> z2_form_factors{};
        std::array<double, max_points> x_prefix{};
        std::array<double, max_points> y_prefix{};
        std::array<double, max_points> z_prefix{};
        std::array<double, max_points * max_points> yz_mix{};
        std::array<double, max_points * max_points> xz_mix{};
        std::array<double, max_points * max_points> xy_mix{};

        const auto mix_pair = [](double a1, double a2, double b1, double b2) {
          return (a1 * b1 + a2 * b2) / 3.0 + (a2 * b1 + a1 * b2) / 6.0;
        };

        const auto fill_support = [](double current_position, int current_cell,
                                     double previous_position,
                                     int previous_cell, int begin, int count,
                                     auto &current_weights,
                                     auto &previous_weights,
                                     auto &prefix_weights) {
          auto fill_weights = [begin, count](double position, int cell,
                                             auto &weights) {
            std::fill(weights.begin(), weights.end(), 0.0);

            const double fraction = position - cell;
            const double one_minus_fraction = 1.0 - fraction;
            const double fraction_sq = fraction * fraction;
            const double fraction_cu = fraction_sq * fraction;
            const double one_minus_fraction_sq =
                one_minus_fraction * one_minus_fraction;
            const double one_minus_fraction_cu =
                one_minus_fraction_sq * one_minus_fraction;

            const std::array<double, 4> local_weights{
                one_minus_fraction_cu / 6.0,
                (3.0 * fraction_cu - 6.0 * fraction_sq + 4.0) / 6.0,
                (-3.0 * fraction_cu + 3.0 * fraction_sq + 3.0 * fraction +
                 1.0) /
                    6.0,
                fraction_cu / 6.0,
            };

            const int support_begin = cell - (particle_radius - 1);
            for (int local_offset = 0; local_offset < 4; ++local_offset) {
              const int union_offset = support_begin - begin + local_offset;
              if (union_offset >= 0 && union_offset < count)
                weights[union_offset] = local_weights[local_offset];
            }
          };

          fill_weights(current_position, current_cell, current_weights);
          fill_weights(previous_position, previous_cell, previous_weights);

          double cumulative = 0.0;
          for (int offset = 0; offset < count; ++offset) {
            cumulative += current_weights[offset] - previous_weights[offset];
            prefix_weights[offset] = cumulative;
          }
        };

        auto *const jx = J.x.v.data();
        auto *const jy = J.y.v.data();
        auto *const jz = J.z.v.data();

        constexpr int jx_sx = Layout<NX, NY, NZ, Component::Jx>::sx;
        constexpr int jx_sy = Layout<NX, NY, NZ, Component::Jx>::sy;
        constexpr int jy_sx = Layout<NX, NY, NZ, Component::Jy>::sx;
        constexpr int jy_sy = Layout<NX, NY, NZ, Component::Jy>::sy;
        constexpr int jz_sx = Layout<NX, NY, NZ, Component::Jz>::sx;
        constexpr int jz_sy = Layout<NX, NY, NZ, Component::Jz>::sy;

        const int start{workers[thread_id].start};
        const int end{workers[thread_id].end};

        for (int idx{start}; idx < end; ++idx) {
          const Particle p{particles[idx]};

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

          const int is = std::max(0, std::min(i, pi) - (particle_radius - 1));
          const int js = std::max(0, std::min(j, pj) - (particle_radius - 1));
          const int ks = std::max(0, std::min(k, pk) - (particle_radius - 1));

          const int ie =
              std::min(J.x.nx() - 1, std::max(i, pi) + particle_radius);
          const int je =
              std::min(J.y.ny() - 1, std::max(j, pj) + particle_radius);
          const int ke =
              std::min(J.z.nz() - 1, std::max(k, pk) + particle_radius);

          const int x_count = ie - is + 1;
          const int y_count = je - js + 1;
          const int z_count = ke - ks + 1;

          fill_support(px, i, ppx, pi, is, x_count, x1_form_factors,
                       x2_form_factors, x_prefix);
          fill_support(py, j, ppy, pj, js, y_count, y1_form_factors,
                       y2_form_factors, y_prefix);
          fill_support(pz, k, ppz, pk, ks, z_count, z1_form_factors,
                       z2_form_factors, z_prefix);

          for (int y_offset = 0; y_offset < y_count; ++y_offset)
            for (int z_offset = 0; z_offset < z_count; ++z_offset)
              yz_mix[y_offset * max_points + z_offset] = mix_pair(
                  y1_form_factors[y_offset], y2_form_factors[y_offset],
                  z1_form_factors[z_offset], z2_form_factors[z_offset]);

          for (int x_offset = 0; x_offset < x_count; ++x_offset)
            for (int z_offset = 0; z_offset < z_count; ++z_offset)
              xz_mix[x_offset * max_points + z_offset] = mix_pair(
                  x1_form_factors[x_offset], x2_form_factors[x_offset],
                  z1_form_factors[z_offset], z2_form_factors[z_offset]);

          for (int x_offset = 0; x_offset < x_count; ++x_offset)
            for (int y_offset = 0; y_offset < y_count; ++y_offset)
              xy_mix[x_offset * max_points + y_offset] = mix_pair(
                  x1_form_factors[x_offset], x2_form_factors[x_offset],
                  y1_form_factors[y_offset], y2_form_factors[y_offset]);

          const double move_co = -p.q / Config::dt;

          for (int x_offset = 0; x_offset < x_count; ++x_offset) {
            const double x_scale = move_co * x_prefix[x_offset];
            const int x_base = (is + x_offset) * jx_sx;

            for (int y_offset = 0; y_offset < y_count; ++y_offset) {
              const int yz_row = y_offset * max_points;
              const int row_base = x_base + (js + y_offset) * jx_sy + ks;

              for (int z_offset = 0; z_offset < z_count; ++z_offset)
                jx[row_base + z_offset] += x_scale * yz_mix[yz_row + z_offset];
            }
          }

          for (int x_offset = 0; x_offset < x_count; ++x_offset) {
            const int xz_row = x_offset * max_points;
            const int x_base = (is + x_offset) * jy_sx;

            for (int y_offset = 0; y_offset < y_count; ++y_offset) {
              const double y_scale = move_co * y_prefix[y_offset];
              const int row_base = x_base + (js + y_offset) * jy_sy + ks;

              for (int z_offset = 0; z_offset < z_count; ++z_offset)
                jy[row_base + z_offset] += y_scale * xz_mix[xz_row + z_offset];
            }
          }

          for (int x_offset = 0; x_offset < x_count; ++x_offset) {
            const int x_base = (is + x_offset) * jz_sx;
            const int xy_row = x_offset * max_points;

            for (int y_offset = 0; y_offset < y_count; ++y_offset) {
              const double xy_scale = xy_mix[xy_row + y_offset];
              const int row_base = x_base + (js + y_offset) * jz_sy + ks;

              for (int z_offset = 0; z_offset < z_count; ++z_offset)
                jz[row_base + z_offset] +=
                    move_co * z_prefix[z_offset] * xy_scale;
            }
          }
        }
      }

      {
        std::lock_guard<std::mutex> lock(m);
        if (thread_id < active_workers) {
          finished++;
          if (finished == active_workers)
            done_cv.notify_one();
        }
      }
    }
  }

  std::vector<Worker> workers;
  std::mutex m;
  std::condition_variable state_cv;
  std::condition_variable done_cv;
  bool stop{false};

  int finished{};
  int active_workers{};
  int job_id{};
  Particle *particles;

  ParallelCurrentDepositor(int thread_count = std::thread::hardware_concurrency()) {
    if (thread_count == 0)
      thread_count = 1;
    workers.reserve(thread_count);
    for (int thread_id{0}; thread_id < thread_count; ++thread_id) {
      workers.push_back(Worker{});
      workers.back().thread =
          std::thread([this, thread_id] { worker_loop(thread_id); });
    }
  }

  ~ParallelCurrentDepositor() {
    {
      std::lock_guard<std::mutex> lock(m);
      stop = true;
    }

    state_cv.notify_all();

    for (auto &worker : workers)
      if (worker.thread.joinable())
        worker.thread.join();
  }

  ParallelCurrentDepositor(const ParallelCurrentDepositor &other) = delete;
  ParallelCurrentDepositor &operator=(const ParallelCurrentDepositor &other) = delete;

  void deposit_currents(std::vector<Particle> &particles, JField<NX, NY, NZ> *j) {
    std::fill(j->x.v.begin(), j->x.v.end(), 0);
    std::fill(j->y.v.begin(), j->y.v.end(), 0);
    std::fill(j->z.v.begin(), j->z.v.end(), 0);

    const int max_threads = workers.size();
    const int num_particles = static_cast<int>(particles.size());
    const int active = std::min(num_particles, max_threads);
    
    const int block = (num_particles + active - 1) / active;

    {
      std::lock_guard<std::mutex> lock(m);

      active_workers = active;
      this->particles = particles.data();

      for (int thread_id{}; thread_id < active; ++thread_id) {
        workers[thread_id].start = thread_id * block;
        workers[thread_id].end = std::min(num_particles, (thread_id + 1) * block);
      }

      finished = 0;
      job_id++;
    }

    state_cv.notify_all();

    {
      std::unique_lock<std::mutex> lock(m);
      done_cv.wait(lock, [&] { return finished == active_workers; });
    }

    for (int thread_id{}; thread_id < active; ++thread_id) {
      auto &worker_j = workers[thread_id].J;
      std::transform(j->x.v.begin(), j->x.v.end(), worker_j.x.v.begin(), j->x.v.begin(), std::plus<double>());
      std::transform(j->y.v.begin(), j->y.v.end(), worker_j.y.v.begin(), j->y.v.begin(), std::plus<double>());
      std::transform(j->z.v.begin(), j->z.v.end(), worker_j.z.v.begin(), j->z.v.begin(), std::plus<double>());
    }
  }
};
