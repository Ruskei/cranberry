#pragma once

template <class T>
struct Vec3 {
  T x{}, y{}, z{};

  Vec3<T> operator-(const Vec3<T> &a) {
    return Vec3<T>{
      x - a.x,
      y - a.y,
      z - a.z,
    };
  };
};
using Vec3d = Vec3<double>;
