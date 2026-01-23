#pragma "once"

#include <vtkImageData.h>
#include <vtkImageMagnitude.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkThresholdPoints.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPolyDataWriter.h>

#include "grid.hpp"

template <int N> struct Writer2D {
  std::string name;
  vtkNew<vtkImageData> image;
  int *dims;
  int slice;

  Writer2D(std::string n, int s) : name{n}, slice{s} {
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
  }

  void write_timestep(Grid<N> &grid, double time) const {
    for (auto y{1}; y < dims[1] - 1; ++y)
      for (auto x{1}; x < dims[0] - 1; ++x) {
        float *pixel = static_cast<float *>(image->GetScalarPointer(x, y, 0));
        pixel[0] = grid.ex(x, y, slice);
        pixel[1] = grid.ey(x, y, slice);
        pixel[2] = grid.ez(x, y, slice);
      }

    image->GetPointData()->SetActiveVectors("ImageScalars");

    image->Modified();

    vtkNew<vtkXMLImageDataWriter> writer;
    const std::string file_name =
        "out/" + name + "/efield" + std::to_string(time) + ".vti";
    writer->SetFileName(file_name.c_str());
    writer->SetInputData(image);
    writer->Write();
  }

  void write_all(double max_time, double dt) const {
    std::ofstream output;
    output.open(("out/" + name + ".pvd").c_str());

    output << "<?xml version=\"1.0\"?>" << std::endl;
    output << "<VTKFile type=\"Collection\" version=\"0.1\" "
              "byte_order=\"LittleEndian\">"
           << std::endl;
    output << "  <Collection>" << std::endl;

    for (auto t{0.0}; t < max_time; t += dt) {
      output << "    <DataSet timestep=\"" << std::to_string(t)
             << "\" group=\"\" part=\"0\" file=\"" + name + "/efield"
             << std::to_string(t) << ".vti\" name=\"efield\"/>" << std::endl;
    }

    output << "  </Collection>" << std::endl;
    output << "</VTKFile>" << std::endl;

    output.close();
  }
};

template <int N> struct Writer3D {
  vtkNew<vtkImageData> image;
  int *dims;

  Writer3D() {
    image->SetDimensions(N, N, N);
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
    mkdir("out/3d_animation", 0755);
  }

  void write_timestep(Grid<N> &grid, double time) const {
    for (auto z{1}; z < dims[2] - 1; ++z)
      for (auto y{1}; y < dims[1] - 1; ++y)
        for (auto x{1}; x < dims[0] - 1; ++x) {
          float *pixel = static_cast<float *>(image->GetScalarPointer(x, y, z));
          pixel[0] = grid.ex(x, y, z);
          pixel[1] = grid.ey(x, y, z);
          pixel[2] = grid.ez(x, y, z);
        }

    image->GetPointData()->SetActiveVectors("ImageScalars");

    image->Modified();

    vtkNew<vtkXMLImageDataWriter> writer;
    const std::string file_name =
        "out/3d_animation/efield" + std::to_string(time) + ".vti";
    writer->SetFileName(file_name.c_str());
    writer->SetInputData(image);
    writer->Write();
  }

  void write_all(double max_time, double dt) const {
    std::ofstream output;
    output.open("out/3d_animation.pvd");

    output << "<?xml version=\"1.0\"?>" << std::endl;
    output << "<VTKFile type=\"Collection\" version=\"0.1\" "
              "byte_order=\"LittleEndian\">"
           << std::endl;
    output << "  <Collection>" << std::endl;

    for (auto t{0.0}; t < max_time; t += dt) {
      output << "    <DataSet timestep=\"" << std::to_string(t)
             << "\" group=\"\" part=\"0\" file=\"3d_animation/efield"
             << std::to_string(t) << ".vti\" name=\"efield\"/>" << std::endl;
    }

    output << "  </Collection>" << std::endl;
    output << "</VTKFile>" << std::endl;

    output.close();
  }
};
