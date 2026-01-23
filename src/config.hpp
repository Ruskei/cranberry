#pragma once

#include <cmath>

namespace Config {
inline const double courants = 1.0 / std::sqrt(3.0);
inline const double abcco = (courants - 1.0) / (courants + 1.0);
inline constexpr int size = 40;
inline constexpr int max_time = 300;
inline constexpr double imp0 = 377.0;
inline constexpr double f0 = 0.08;
inline constexpr double t0 = 30.0;
}; // namespace Config
