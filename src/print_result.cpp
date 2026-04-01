#include <iostream>

#include "print_result.hpp"
#include "config.hpp"

inline void print_progress(double frac, int width) {
  const int filled = static_cast<int>(frac * width);

  std::cout << "\r[";
  for (auto i{0}; i < width; ++i)
    std::cout << (i < filled ? '=' : ' ');

  std::cout << "] " << static_cast<int>(frac * 100) << "%";
  std::cout.flush();
}

inline void print_sim_finished() {
  std::cout << "\r\x1b[2K";
  std::cout << std::flush;
  std::cout << "Finished simulation!" << std::endl;
}

void PrintResult::step(double time) {
  if (static_cast<int>(time / Config::dt) % Config::print_interval == 0)
    print_progress(time / Config::max_time, 20);
}
