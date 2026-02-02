#pragma once

#include <algorithm>
#include <cassert>
#include <chrono>
#include <condition_variable>
#include <iostream>
#include <mutex>
#include <ratio>
#include <thread>

#include "runtime_field.hpp"

// NOT safe to call from multiple threads
struct MultigridJacobi {
  // processes [x0, x1)

  struct Worker {
    std::thread thread;
    int x0{}, x1{};
  };

  void worker_loop(unsigned thread_id) {
    int local_job_id{};
    for (;;) {
      {
        std::unique_lock<std::mutex> lock(m);
        state_cv.wait(lock, [&] { return stop || local_job_id != job_id; });
        if (stop)
          return;

        local_job_id = job_id;
      }

      if (!job_out || !job_guess || !job_target)
        return;

      const int sy{job_guess->sy}, sz{job_guess->sz};
      const int yz{sy * sz};

      const int x0 = workers[thread_id].x0;
      const int x1 = workers[thread_id].x1;

      if (x1 > x0) {
        for (auto x{x0}; x < x1; ++x) {
          const int x_idx{x * yz};

          for (auto y{1}; y < sy - 1; ++y) {
            const int y_idx{x_idx + y * sz};

            for (auto z{1}; z < sz - 1; ++z) {
              const int idx = y_idx + z;

              const double jacobi =
                  (job_guess->v[idx - yz] + job_guess->v[idx + yz] +
                   job_guess->v[idx - sz] + job_guess->v[idx + sz] +
                   job_guess->v[idx - 1] + job_guess->v[idx + 1] +
                   job_h2 * job_target->v[idx]) /
                  6.0;

              job_out->v[idx] =
                  (1.0 - job_omega) * job_guess->v[idx] + job_omega * jacobi;
            }
          }
        }
      }

      {
        std::lock_guard<std::mutex> lock(m);
        if (x1 > x0) {
          finished++;
          if (finished == active_workers)
            done_cv.notify_one();
        }
      }
    }
  }

  RuntimeField *job_guess{nullptr};
  RuntimeField *job_out{nullptr};
  const RuntimeField *job_target{nullptr};
  double job_h2{}, job_omega{};

  std::mutex m;
  std::condition_variable state_cv;
  std::condition_variable done_cv;
  bool stop{false};
  int job_id{};

  int active_workers{};
  int finished = 0;

  std::vector<Worker> workers;

  MultigridJacobi(size_t thread_count = std::thread::hardware_concurrency()) {
    if (thread_count == 0)
      thread_count = 1;
    workers.reserve(thread_count);
    for (size_t thread_id{0}; thread_id < thread_count; ++thread_id) {
      workers.push_back(Worker{});
      workers.back().thread =
          std::thread([this, thread_id] { worker_loop(thread_id); });
    }
  }

  ~MultigridJacobi() {
    {
      std::lock_guard<std::mutex> lock(m);
      stop = true;
    }

    state_cv.notify_all();

    for (auto &worker : workers)
      if (worker.thread.joinable())
        worker.thread.join();
  }

  MultigridJacobi(const MultigridJacobi &other) = delete;
  MultigridJacobi &operator=(const MultigridJacobi &other) = delete;

  void enforce_boundary_conditions(RuntimeField &field) {
    const int sx = field.sx, sy = field.sy, sz = field.sz;
    const int yz = sy * sz;

    std::fill_n(field.v.begin(), (sy - 1) * sz + (sz - 1), 0.0);
    std::fill_n(field.v.begin() + (sz - 1) * yz, (sy - 1) * sz + (sz - 1), 0.0);

    for (int x = 0; x < sx; ++x) {
      std::fill_n(field.v.begin() + (x * yz), sz - 1, 0.0);
      std::fill_n(field.v.begin() + (x * yz) + (sy - 1) * sz, sz - 1, 0.0);
    }

    for (int x = 0; x < sx; ++x)
      for (int y = 0; y < sy; ++y)
        field.v[x * yz + y * sz] = field.v[x * yz + y * sz + (sz - 1)] = 0.0;
  }

  void smooth_serial(RuntimeField &guess, const RuntimeField &target,
                     double scale_factor, double omega = 2.0 / 3.0) {
    const int sx = guess.sx, sy = guess.sy, sz = guess.sz;
    const int yz = sy * sz;
    const double h2 = scale_factor * scale_factor;

    RuntimeField new_guess{sx, sy, sz};
    for (int x = 1; x < sx - 1; ++x) {
      const int xb = x * yz;
      for (int y = 1; y < sy - 1; ++y) {
        const int yb = xb + y * sz;
        for (int z = 1; z < sz - 1; ++z) {
          const int i = yb + z;

          const double jacobi =
              (guess.v[i - yz] + guess.v[i + yz] + guess.v[i - sz] +
               guess.v[i + sz] + guess.v[i - 1] + guess.v[i + 1] +
               h2 * target.v[i]) /
              6.0;

          new_guess.v[i] = (1.0 - omega) * guess.v[i] + omega * jacobi;
        }
      }
    }

    guess.v.swap(new_guess.v);
  }

  void smooth(RuntimeField &guess, const RuntimeField &target,
              RuntimeField &temp, double scale_factor,
              double omega = 2.0 / 3.0) {
    auto start = std::chrono::high_resolution_clock::now();

    assert(guess.sx == target.sx && guess.sy == target.sy &&
           guess.sz == target.sz);
    const int sx = guess.sx, sy = guess.sy, sz = guess.sz;

    const long long volume = 1LL * sx * sy * sz;
    const long long threshold = 1LL * 64 * 64 * 64;

    if (volume <= threshold || workers.size() == 1) {
      smooth_serial(guess, target, scale_factor, omega);
      return;
    }

    const double h2 = scale_factor * scale_factor;

    const int x0 = 1;
    const int x1 = sx - 1;
    const int x_range = x1 - x0;
    if (x_range <= 0)
      return;

    const int max_threads = static_cast<int>(workers.size());
    const int active = std::min(x_range, static_cast<int>(max_threads));

    const int block =
        (x_range + static_cast<int>(active) - 1) / static_cast<int>(active);

    {
      std::lock_guard<std::mutex> lock(m);

      job_out = &temp;
      job_guess = &guess;
      job_target = &target;
      job_h2 = h2;
      job_omega = omega;

      active_workers = active;

      for (int t{0}; t < active; ++t) {
        const int local_x0 = x0 + t * block;
        const int local_x1 = std::min(x1, local_x0 + block);
        workers[t].x0 = local_x0;
        workers[t].x1 = local_x1;
      }

      for (int t{active}; t < max_threads; ++t) {
        workers[t].x0 = 0;
        workers[t].x1 = 0;
      }

      finished = 0;
      job_id++;
    }

    auto setup = std::chrono::high_resolution_clock::now();
    state_cv.notify_all();

    {
      std::unique_lock<std::mutex> lock(m);
      done_cv.wait(lock, [&] { return finished == active_workers; });
    }

    auto solving = std::chrono::high_resolution_clock::now();

    enforce_boundary_conditions(temp);
    guess.v.swap(temp.v);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> setup_duration{setup - start};
    std::chrono::duration<double, std::milli> solving_duration{solving - setup};
    std::chrono::duration<double, std::milli> saving_duration{finish - solving};
    std::chrono::duration<double, std::milli> total_duration{finish - start};
    // std::cout << "took " << total_duration.count() << " ms\n"
    //           << "  setup= " << setup_duration.count() << "\n"
    //           << "  solving= " << solving_duration.count() << "\n"
    //           << "  saving=" << saving_duration.count() << std::endl;
  }
};
