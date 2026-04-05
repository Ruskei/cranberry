#pragma once

#include <condition_variable>
#include <mutex>
#include <thread>

#include "fdtd_types.hpp"
#include "config.hpp"
#include "particles.hpp"

template <int NX, int NY, int NZ> struct ParallelParticlePusher {
  struct Worker {
    std::thread thread;
    int start{}, end{};
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
        const double epsilon = 1e-12;
        const double dt = Config::dt;

        const int start{workers[thread_id].start};
        const int end{workers[thread_id].end};

        for (int idx{start}; idx < end; ++idx) {
          Particle *p = &particles[idx];
          const double nx = p->px;
          const double ny = p->py;
          const double nz = p->pz;

          double vx = p->vx;
          double vy = p->vy;
          double vz = p->vz;

          double vv = vx * vx + vy * vy + vz * vz;
          if (!std::isfinite(vv) || vv >= 1.0 - epsilon) {
            double scale = std::sqrt((1.0 - 1e-12) / std::max(vv, epsilon));
            vx *= scale;
            vy *= scale;
            vz *= scale;
            vv = vx * vx + vy * vy + vz * vz;
          }

          const double q = p->q;
          const double m = p->m;

          double lorentz = 1.0 / std::sqrt(1.0 - vv);

          const double exv = E->x.interpolate_at(nx, ny, nz);
          const double eyv = E->y.interpolate_at(nx, ny, nz);
          const double ezv = E->z.interpolate_at(nx, ny, nz);

          const double hxv = H->x.interpolate_at(nx, ny, nz);
          const double hyv = H->y.interpolate_at(nx, ny, nz);
          const double hzv = H->z.interpolate_at(nx, ny, nz);

          const double ux = lorentz * vx;
          const double uy = lorentz * vy;
          const double uz = lorentz * vz;

          const double u_minus_x = ux + (q * dt / (2.0 * m)) * exv;
          const double u_minus_y = uy + (q * dt / (2.0 * m)) * eyv;
          const double u_minus_z = uz + (q * dt / (2.0 * m)) * ezv;
          const double umum = (u_minus_x * u_minus_x + u_minus_y * u_minus_y +
                               u_minus_z * u_minus_z);

          const double lorentz_bar = 1.0 / sqrt(1.0 + umum);
          const double tx = q * dt / (2.0 * m) * hxv * lorentz_bar;
          const double ty = q * dt / (2.0 * m) * hyv * lorentz_bar;
          const double tz = q * dt / (2.0 * m) * hzv * lorentz_bar;
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

          const double u_next_x = u_plus_x + q * dt / (2 * m) * exv;
          const double u_next_y = u_plus_y + q * dt / (2 * m) * eyv;
          const double u_next_z = u_plus_z + q * dt / (2 * m) * ezv;
          const double unun =
              (u_next_x * u_next_x + u_next_y * u_next_y + u_next_z * u_next_z);

          const double lorentz_next = std::sqrt(1.0 + unun);
          p->vx = u_next_x / lorentz_next;
          p->vy = u_next_y / lorentz_next;
          p->vz = u_next_z / lorentz_next;

          p->p_prev_x = p->px;
          p->p_prev_y = p->py;
          p->p_prev_z = p->pz;

          p->px += p->vx * dt;
          p->py += p->vy * dt;
          p->pz += p->vz * dt;
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
  EField<NX, NY, NZ> *E{nullptr};
  HField<NX, NY, NZ> *H{nullptr};

  ParallelParticlePusher(int thread_count = std::thread::hardware_concurrency()) {
    if (thread_count == 0)
      thread_count = 1;
    workers.reserve(thread_count);
    for (int thread_id{0}; thread_id < thread_count; ++thread_id) {
      workers.push_back(Worker{});
      workers.back().thread =
          std::thread([this, thread_id] { worker_loop(thread_id); });
    }
  }

  ~ParallelParticlePusher() {
    {
      std::lock_guard<std::mutex> lock(m);
      stop = true;
    }

    state_cv.notify_all();

    for (auto &worker : workers)
      if (worker.thread.joinable())
        worker.thread.join();
  }

  ParallelParticlePusher(const ParallelParticlePusher &other) = delete;
  ParallelParticlePusher &operator=(const ParallelParticlePusher &other) = delete;

  void push_particles(EField<NX, NY, NZ> *e, HField<NX, NY, NZ> *h, std::vector<Particle> &particles) {
    const int max_threads = workers.size();
    const int num_particles = static_cast<int>(particles.size());
    const int active = std::min(num_particles, max_threads);
    
    const int block = (num_particles + active - 1) / active;

    {
      std::lock_guard<std::mutex> lock(m);


      active_workers = active;
      E = e;
      H = h;
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
  }
};
