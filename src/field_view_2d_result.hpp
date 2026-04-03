#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <utility>

#include <vtkImageData.h>
#include <vtkImageMagnitude.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkThresholdPoints.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPolyDataWriter.h>

#include "config.hpp"
#include "fdtd_types.hpp"

enum Axis { X, Y, Z };
enum WhichField { E, J, H };

template <int NX, int NY, int NZ> struct FieldView2D {
  using Shape = GridShape<NX, NY, NZ>;
  std::string name;
  WhichField field;
  Axis axis;
  int slice{};
  int interval{1};
  double threshold{0.995};

  vtkNew<vtkImageData> image;
  int *dims;
  int num_saved{};

  static constexpr int plane_width(Axis axis) {
    switch (axis) {
    case Axis::X:
      return NY - 1;
    case Axis::Y:
      return NX - 1;
    case Axis::Z:
      return NX - 1;
    }

    return 0;
  }

  static constexpr int plane_height(Axis axis) {
    switch (axis) {
    case Axis::X:
      return NZ - 1;
    case Axis::Y:
      return NZ - 1;
    case Axis::Z:
      return NY - 1;
    }

    return 0;
  }

  static constexpr int slice_extent(Axis axis) {
    switch (axis) {
    case Axis::X:
      return NX - 1;
    case Axis::Y:
      return NY - 1;
    case Axis::Z:
      return NZ - 1;
    }

    return 0;
  }

  template <class VectorField>
  std::array<double, 3> sample(const VectorField &field, int a, int b) const {
    if (axis == Axis::X)
      return {field.x(slice, a, b), field.y(slice, a, b), field.z(slice, a, b)};
    if (axis == Axis::Y)
      return {field.x(a, slice, b), field.y(a, slice, b), field.z(a, slice, b)};

    return {field.x(a, b, slice), field.y(a, b, slice), field.z(a, b, slice)};
  }

  FieldView2D(std::string name, WhichField field, Axis axis, int slice,
              int interval)
      : name{std::move(name)}, field{field}, axis{axis}, interval{interval} {
    this->slice = std::clamp(slice, 0, slice_extent(axis) - 1);

    image->SetDimensions(plane_width(axis), plane_height(axis), 1);
    image->AllocateScalars(VTK_FLOAT, 3);

    dims = image->GetDimensions();

    for (auto z{0}; z < dims[2]; ++z)
      for (auto y{0}; y < dims[1]; ++y)
        for (auto x{0}; x < dims[0]; ++x) {
          float *pixel = static_cast<float *>(image->GetScalarPointer(x, y, z));
          pixel[0] = 0;
          pixel[1] = 0;
          pixel[2] = 0;
        }

    mkdir("out", 0755);
    mkdir(("out/" + name).c_str(), 0755);
    mkdir(("out/" + name + "/field").c_str(), 0755);
  }

  void step(
      double time,
      const EField<NX, NY, NZ> &E,
      const JField<NX, NY, NZ> &J,
      const HField<NX, NY, NZ> &H
    ) {
    if (static_cast<int>(std::round(time / Config::dt)) % interval != 0) return;

    std::vector<double> mags;
    mags.reserve((dims[1] - 2) * (dims[0] - 2));

    for (auto y{1}; y < dims[1] - 1; ++y)
      for (auto x{1}; x < dims[0] - 1; ++x) {
        const std::array<double, 3> values =
            (field == WhichField::E)
                ? sample(E, x, y)
                : (field == WhichField::J) ? sample(J, x, y) : sample(H, x, y);
        const double vx = values[0];
        const double vy = values[1];
        const double vz = values[2];
        mags.push_back(std::sqrt(vx * vx + vy * vy + vz * vz));
      }

    const int k = std::floor(threshold * (mags.size() - 1));
    std::nth_element(mags.begin(), mags.begin() + k, mags.end());
    const double cutoff = mags[k];
    for (auto y{1}; y < dims[1] - 1; ++y)
      for (auto x{1}; x < dims[0] - 1; ++x) {
        float *pixel = static_cast<float *>(image->GetScalarPointer(x, y, 0));
        std::array<double, 3> values =
            (field == WhichField::E)
                ? sample(E, x, y)
                : (field == WhichField::J) ? sample(J, x, y) : sample(H, x, y);
        double vx = values[0];
        double vy = values[1];
        double vz = values[2];

        const double m = std::sqrt(vx * vx + vy * vy + vz * vz);
        if (m > cutoff) {
          vx /= (m / cutoff);
          vy /= (m / cutoff);
          vz /= (m / cutoff);
        }

        pixel[0] = vx;
        pixel[1] = vy;
        pixel[2] = vz;
      }

    image->GetPointData()->SetActiveVectors("ImageScalars");

    image->Modified();

    vtkNew<vtkXMLImageDataWriter> writer;
    const std::string file_name =
        "out/" + name + "/field/" + std::to_string(num_saved++) + ".vti";
    writer->SetFileName(file_name.c_str());
    writer->SetInputData(image);
    writer->Write();
  }

  void finish() {
    std::ofstream output;
    output.open(("out/" + name + "/recording.pvd").c_str());

    output << "<?xml version=\"1.0\"?>" << std::endl;
    output << "<VTKFile type=\"Collection\" version=\"0.1\" "
              "byte_order=\"LittleEndian\">"
           << std::endl;
    output << "  <Collection>" << std::endl;

    for (auto i{0}; i < num_saved; i++) {
      output << "    <DataSet timestep=\"" << std::to_string(i)
             << "\" group=\"\" part=\"0\" file=\"field/"
             << std::to_string(i) << ".vti\" name=\"field_view_2d\"/>" << std::endl;
    }

    output << "  </Collection>" << std::endl;
    output << "</VTKFile>" << std::endl;

    output.close();
  }
};
