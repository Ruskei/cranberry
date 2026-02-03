#include <algorithm>
#include <condition_variable>
#include <mutex>
#include <thread>

#include "fdtd_types.hpp"

template <int N> struct ParallelFDTD {
  enum Type { none, e, h };

  struct Worker {
    std::thread thread;
    int x_x0{}, x_x1{}, y_x0{}, y_x1{}, z_x0{}, z_x1{};
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

      if (type == Type::h) {
        if (!E || !H)
          return;

        const int x_x0{workers[thread_id].x_x0};
        const int x_x1{workers[thread_id].x_x1};

        const int y_x0{workers[thread_id].y_x0};
        const int y_x1{workers[thread_id].y_x1};

        const int z_x0{workers[thread_id].z_x0};
        const int z_x1{workers[thread_id].z_x1};

        if (thread_id < active_workers) {
          for (auto x{x_x0}; x < x_x1; ++x)
            for (auto y{0}; y < H->x.ny(); ++y)
              for (auto z{0}; z < H->x.nz(); ++z) {
                const double dey_dz = E->y(x, y, z + 1) - E->y(x, y, z);
                const double dez_dy = E->z(x, y + 1, z) - E->z(x, y, z);
                H->x(x, y, z) = H->cxh(x, y, z) * H->x(x, y, z) +
                                0.5 * H->cxe(x, y, z) * (dey_dz - dez_dy);
              }

          for (auto x{y_x0}; x < y_x1; ++x)
            for (auto y{0}; y < H->y.ny(); ++y)
              for (auto z{0}; z < H->y.nz(); ++z) {
                const double dex_dz = E->x(x, y, z + 1) - E->x(x, y, z);
                const double dez_dx = E->z(x + 1, y, z) - E->z(x, y, z);
                H->y(x, y, z) = H->cyh(x, y, z) * H->y(x, y, z) +
                                0.5 * H->cye(x, y, z) * (dez_dx - dex_dz);
              }

          for (auto x{z_x0}; x < z_x1; ++x)
            for (auto y{0}; y < H->z.ny(); ++y)
              for (auto z{0}; z < H->z.nz(); ++z) {
                const double dex_dy = E->x(x, y + 1, z) - E->x(x, y, z);
                const double dey_dx = E->y(x + 1, y, z) - E->y(x, y, z);
                H->z(x, y, z) = H->czh(x, y, z) * H->z(x, y, z) +
                                0.5 * H->cze(x, y, z) * (dex_dy - dey_dx);
              }
        }
      } else if (type == Type::e) {
        if (!E || !H || !J)
          return;

        const int x_x0{workers[thread_id].x_x0};
        const int x_x1{workers[thread_id].x_x1};

        const int y_x0{workers[thread_id].y_x0};
        const int y_x1{workers[thread_id].y_x1};

        const int z_x0{workers[thread_id].z_x0};
        const int z_x1{workers[thread_id].z_x1};

        if (thread_id < active_workers) {
          for (auto x{x_x0}; x < x_x1; ++x)
            for (auto y{1}; y < E->x.ny() - 1; ++y)
              for (auto z{1}; z < E->x.nz() - 1; ++z) {
                const double dhz_dy = H->z(x, y, z) - H->z(x, y - 1, z);
                const double dhy_dz = H->y(x, y, z) - H->y(x, y, z - 1);
                const double j = J->x(x, y, z);
                E->x(x, y, z) = E->cxe(x, y, z) * E->x(x, y, z) +
                               0.5 * E->cxh(x, y, z) * (dhz_dy - dhy_dz - j);
              }

          for (auto x{y_x0}; x < y_x1; ++x)
            for (auto y{0}; y < E->y.ny(); ++y)
              for (auto z{1}; z < E->y.nz() - 1; ++z) {
                const double dhx_dhz = H->x(x, y, z) - H->x(x, y, z - 1);
                const double dhz_dhx = H->z(x, y, z) - H->z(x - 1, y, z);
                const double j = J->y(x, y, z);
                E->y(x, y, z) = E->cye(x, y, z) * E->y(x, y, z) +
                               0.5 * E->cyh(x, y, z) * (dhx_dhz - dhz_dhx - j);
              }

          for (auto x{z_x0}; x < z_x1; ++x)
            for (auto y{1}; y < E->z.ny() - 1; ++y)
              for (auto z{0}; z < E->z.nz(); ++z) {
                const double dhy_dhx = H->y(x, y, z) - H->y(x - 1, y, z);
                const double dhx_dhy = H->x(x, y, z) - H->x(x, y - 1, z);
                const double j = J->z(x, y, z);
                E->z(x, y, z) = E->cze(x, y, z) * E->z(x, y, z) +
                               0.5 * E->czh(x, y, z) * (dhy_dhx - dhx_dhy - j);
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
  Type type{none};
  int job_id{};
  EField<N> *E{nullptr};
  HField<N> *H{nullptr};
  JField<N> *J{nullptr};

  ParallelFDTD(int thread_count = std::thread::hardware_concurrency()) {
    if (thread_count == 0)
      thread_count = 1;
    workers.reserve(thread_count);
    for (int thread_id{0}; thread_id < thread_count; ++thread_id) {
      workers.push_back(Worker{});
      workers.back().thread =
          std::thread([this, thread_id] { worker_loop(thread_id); });
    }
  }

  ~ParallelFDTD() {
    {
      std::lock_guard<std::mutex> lock(m);
      stop = true;
    }

    state_cv.notify_all();

    for (auto &worker : workers)
      if (worker.thread.joinable())
        worker.thread.join();
  }

  ParallelFDTD(const ParallelFDTD &other) = delete;
  ParallelFDTD &operator=(const ParallelFDTD &other) = delete;

  void half_update_h_serial(const EField<N> &E, HField<N> &H) {
    for (auto x{0}; x < H.x.nx(); ++x)
      for (auto y{0}; y < H.x.ny(); ++y)
        for (auto z{0}; z < H.x.nz(); ++z) {
          const double dey_dz = E.y(x, y, z + 1) - E.y(x, y, z);
          const double dez_dy = E.z(x, y + 1, z) - E.z(x, y, z);
          H.x(x, y, z) = H.cxh(x, y, z) * H.x(x, y, z) +
                         0.5 * H.cxe(x, y, z) * (dey_dz - dez_dy);
        }

    for (auto x{0}; x < H.y.nx(); ++x)
      for (auto y{0}; y < H.y.ny(); ++y)
        for (auto z{0}; z < H.y.nz(); ++z) {
          const double dex_dz = E.x(x, y, z + 1) - E.x(x, y, z);
          const double dez_dx = E.z(x + 1, y, z) - E.z(x, y, z);
          H.y(x, y, z) = H.cyh(x, y, z) * H.y(x, y, z) +
                         0.5 * H.cye(x, y, z) * (dez_dx - dex_dz);
        }

    for (auto x{0}; x < H.z.nx(); ++x)
      for (auto y{0}; y < H.z.ny(); ++y)
        for (auto z{0}; z < H.z.nz(); ++z) {
          const double dex_dy = E.x(x, y + 1, z) - E.x(x, y, z);
          const double dey_dx = E.y(x + 1, y, z) - E.y(x, y, z);
          H.z(x, y, z) = H.czh(x, y, z) * H.z(x, y, z) +
                         0.5 * H.cze(x, y, z) * (dex_dy - dey_dx);
        }
  }

  void half_update_e_serial(EField<N> &E, const HField<N> &H,
                            const JField<N> &J) {
    for (auto x{0}; x < E.x.nx(); ++x)
      for (auto y{1}; y < E.x.ny() - 1; ++y)
        for (auto z{1}; z < E.x.nz() - 1; ++z) {
          const double dhz_dy = H.z(x, y, z) - H.z(x, y - 1, z);
          const double dhy_dz = H.y(x, y, z) - H.y(x, y, z - 1);
          const double j = J.x(x, y, z);
          E.x(x, y, z) = E.cxe(x, y, z) * E.x(x, y, z) +
                         0.5 * E.cxh(x, y, z) * (dhz_dy - dhy_dz - j);
        }

    for (auto x{1}; x < E.y.nx() - 1; ++x)
      for (auto y{0}; y < E.y.ny(); ++y)
        for (auto z{1}; z < E.y.nz() - 1; ++z) {
          const double dhx_dhz = H.x(x, y, z) - H.x(x, y, z - 1);
          const double dhz_dhx = H.z(x, y, z) - H.z(x - 1, y, z);
          const double j = J.y(x, y, z);
          E.y(x, y, z) = E.cye(x, y, z) * E.y(x, y, z) +
                         0.5 * E.cyh(x, y, z) * (dhx_dhz - dhz_dhx - j);
        }

    for (auto x{1}; x < E.z.nx() - 1; ++x)
      for (auto y{1}; y < E.z.ny() - 1; ++y)
        for (auto z{0}; z < E.z.nz(); ++z) {
          const double dhy_dhx = H.y(x, y, z) - H.y(x - 1, y, z);
          const double dhx_dhy = H.x(x, y, z) - H.x(x, y - 1, z);
          const double j = J.z(x, y, z);
          E.z(x, y, z) = E.cze(x, y, z) * E.z(x, y, z) +
                         0.5 * E.czh(x, y, z) * (dhy_dhx - dhx_dhy - j);
        }
  }

  void half_update(Type t, EField<N> *e, HField<N> *h, JField<N> *j,
                   int x_x_min, int x_x_max, int y_x_min, int y_x_max,
                   int z_x_min, int z_x_max) {
    const int x_x_range = x_x_max - x_x_min;
    const int y_x_range = y_x_max - y_x_min;
    const int z_x_range = z_x_max - z_x_min;

    const int max_threads = workers.size();
    const int active = std::min(x_x_range, max_threads);

    const int x_block = (x_x_range + active - 1) / active;
    const int y_block = (y_x_range + active - 1) / active;
    const int z_block = (z_x_range + active - 1) / active;

    {
      std::lock_guard<std::mutex> lock(m);

      type = t;

      E = e;
      H = h;
      J = j;

      active_workers = active;

      for (int thread_id{}; thread_id < active; ++thread_id) {
        const int local_x_x0 = x_x_min + thread_id * x_block;
        const int local_x_x1 = std::min(x_x_range, local_x_x0 + x_block);
        workers[thread_id].x_x0 = local_x_x0;
        workers[thread_id].x_x1 = local_x_x1;

        const int local_y_x0 = y_x_min + thread_id * y_block;
        const int local_y_x1 = std::min(y_x_range, local_y_x0 + y_block);
        workers[thread_id].y_x0 = local_y_x0;
        workers[thread_id].y_x1 = local_y_x1;

        const int local_z_x0 = z_x_min + thread_id * z_block;
        const int local_z_x1 = std::min(z_x_range, local_z_x0 + z_block);
        workers[thread_id].z_x0 = local_z_x0;
        workers[thread_id].z_x1 = local_z_x1;
      }

      for (int thread_id{active}; thread_id < max_threads; ++thread_id) {
        workers[thread_id].x_x0 = workers[thread_id].x_x1 = 0;
        workers[thread_id].y_x0 = workers[thread_id].y_x1 = 0;
        workers[thread_id].z_x0 = workers[thread_id].z_x1 = 0;
      }

      finished = 0;
      job_id++;
    }

    state_cv.notify_all();

    {
      std::unique_lock<std::mutex> lock(m);
      done_cv.wait(lock, [&] { return finished == active_workers; });
      type = Type::none;
    }
  }

  void half_update_h(EField<N> &e, HField<N> &h) {
    const long long threshold = 1LL * 64 * 64 * 64;
    if (h.x.v.size() < threshold || workers.size() == 1) {
      half_update_h_serial(e, h);
      return;
    }

    half_update(Type::h, &e, &h, nullptr, 0, h.x.nx(), 0, h.y.nx(), 0,
                h.z.nx());
  }

  void half_update_e(EField<N> &e, HField<N> &h, JField<N> &j) {
    const long long threshold = 1LL * 64 * 64 * 64;
    if (h.x.v.size() < threshold || workers.size() == 1) {
      half_update_e_serial(e, h, j);
      return;
    }

    half_update(Type::e, &e, &h, &j, 0, e.x.nx(), 1, e.y.nx() - 1, 1,
                e.z.nx() - 1);
  }
};
