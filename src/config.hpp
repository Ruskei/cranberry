#pragma once

namespace Config { // max courant's is 1 / sqrt(3)
inline constexpr int print_interval = 10;

inline const double dt = 0.5;
inline const double courants = dt;
inline const double abcco = (courants - 1.0) / (courants + 1.0);
inline constexpr int size = 40;
inline constexpr double max_time = 30.0;
inline constexpr double imp0 = 1.0;
inline constexpr double f0 = 0.08;
inline constexpr double t0 = 30.0;
inline constexpr double A = 100.0;
}; // namespace Config
