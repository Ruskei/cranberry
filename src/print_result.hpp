#pragma once

struct PrintResult {
  int print_interval{1};
  void step(double time);
};
