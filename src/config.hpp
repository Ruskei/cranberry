#pragma once

namespace Config { // max courant's is 1 / sqrt(3)
namespace config_detail {
constexpr bool is_power_of_two(int value) {
  return value > 0 && (value & (value - 1)) == 0;
}

constexpr bool is_valid_extent(int value) {
  return value >= 5 && is_power_of_two(value - 1);
}
} // namespace config_detail

inline constexpr int print_interval = 10;

inline const double dt = 0.5;
inline const double courants = dt;
inline const double abcco = (courants - 1.0) / (courants + 1.0);
inline constexpr int nx = 65;
inline constexpr int ny = 129;
inline constexpr int nz = 65;
inline constexpr double max_time = 20.0;
inline constexpr double imp0 = 1.0;
inline constexpr double f0 = 0.08;
inline constexpr double t0 = 30.0;
inline constexpr double A = 100.0;

static_assert(config_detail::is_valid_extent(nx),
              "Config::nx must be 2^n + 1 and >= 5");
static_assert(config_detail::is_valid_extent(ny),
              "Config::ny must be 2^n + 1 and >= 5");
static_assert(config_detail::is_valid_extent(nz),
              "Config::nz must be 2^n + 1 and >= 5");
}; // namespace Config
