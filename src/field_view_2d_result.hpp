#pragma once

#include <string>

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

template <int N> struct FieldView2D {
  std::string name;
  WhichField field;
  Axis axis;
  int slice{};
  int interval{1};
  double threshold{0.995};

  vtkNew<vtkImageData> image;
  int *dims;
  int num_saved{};

  FieldView2D(std::string name, WhichField field, Axis axis, int slice,
              int interval)
      : name{name}, field{field}, axis{axis}, slice{slice}, interval{interval} {
    image->SetDimensions(N, N, 1);
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

  void step(double time, EField<N> &E, JField<N> &J, HField<N> &H) {
    if (static_cast<int>(std::round(time / Config::dt)) % interval != 0) return;

    std::vector<double> mags;
    mags.reserve((dims[1] - 2) * (dims[0] - 2));

    for (auto y{1}; y < dims[1] - 1; ++y)
      for (auto x{1}; x < dims[0] - 1; ++x) {
        double vx, vy, vz;
        if (field == WhichField::E) {
          if (axis == Axis::X) {
            vx = E.x(slice, x, y);
            vy = E.y(slice, x, y);
            vz = E.z(slice, x, y);
          } else if (axis == Axis::Y) {
            vx = E.x(x, slice, y);
            vy = E.y(x, slice, y);
            vz = E.z(x, slice, y);
          } else {
            vx = E.x(x, y, slice);
            vy = E.y(x, y, slice);
            vz = E.z(x, y, slice);
          }
        } else if (field == WhichField::J) {
          if (axis == Axis::X) {
            vx = J.x(slice, x, y);
            vy = J.y(slice, x, y);
            vz = J.z(slice, x, y);
          } else if (axis == Axis::Y) {
            vx = J.x(x, slice, y);
            vy = J.y(x, slice, y);
            vz = J.z(x, slice, y);
          } else {
            vx = J.x(x, y, slice);
            vy = J.y(x, y, slice);
            vz = J.z(x, y, slice);
          }
        } else {
          if (axis == Axis::X) {
            vx = H.x(slice, x, y);
            vy = H.y(slice, x, y);
            vz = H.z(slice, x, y);
          } else if (axis == Axis::Y) {
            vx = H.x(x, slice, y);
            vy = H.y(x, slice, y);
            vz = H.z(x, slice, y);
          } else {
            vx = H.x(x, y, slice);
            vy = H.y(x, y, slice);
            vz = H.z(x, y, slice);
          }
        }
        mags.push_back(std::sqrt(vx * vx + vy * vy + vz * vz));
      }

    const int k = std::floor(threshold * (mags.size() - 1));
    std::nth_element(mags.begin(), mags.begin() + k, mags.end());
    const double cutoff = mags[k];
    for (auto y{1}; y < dims[1] - 1; ++y)
      for (auto x{1}; x < dims[0] - 1; ++x) {
        float *pixel = static_cast<float *>(image->GetScalarPointer(x, y, 0));
        double vx, vy, vz;

        if (field == WhichField::E) {
          if (axis == Axis::X) {
            vx = E.x(slice, x, y);
            vy = E.y(slice, x, y);
            vz = E.z(slice, x, y);
          } else if (axis == Axis::Y) {
            vx = E.x(x, slice, y);
            vy = E.y(x, slice, y);
            vz = E.z(x, slice, y);
          } else {
            vx = E.x(x, y, slice);
            vy = E.y(x, y, slice);
            vz = E.z(x, y, slice);
          }
        } else if (field == WhichField::J) {
          if (axis == Axis::X) {
            vx = J.x(slice, x, y);
            vy = J.y(slice, x, y);
            vz = J.z(slice, x, y);
          } else if (axis == Axis::Y) {
            vx = J.x(x, slice, y);
            vy = J.y(x, slice, y);
            vz = J.z(x, slice, y);
          } else {
            vx = J.x(x, y, slice);
            vy = J.y(x, y, slice);
            vz = J.z(x, y, slice);
          }
        } else {
          if (axis == Axis::X) {
            vx = H.x(slice, x, y);
            vy = H.y(slice, x, y);
            vz = H.z(slice, x, y);
          } else if (axis == Axis::Y) {
            vx = H.x(x, slice, y);
            vy = H.y(x, slice, y);
            vz = H.z(x, slice, y);
          } else {
            vx = H.x(x, y, slice);
            vy = H.y(x, y, slice);
            vz = H.z(x, y, slice);
          }
        }

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
