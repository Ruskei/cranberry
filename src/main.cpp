#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <cmath>
#include <ostream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

#include <vtkImageData.h>
#include <vtkImageMagnitude.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkThresholdPoints.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPolyDataWriter.h>

constexpr int size = 200;

constexpr double imp0 = 377.0;
constexpr int max_time = 100;

int main() {
  std::cout << "Running simulation\n";

  std::array<double, size * max_time> history_ez{};

  std::array<double, size> ez{};
  std::array<double, size> hy{};

  for (auto time{ 0 }; time < max_time; ++time) {
    // update magnetic field
    for (auto mm{ 0 }; mm < size - 1; ++mm)
      hy[mm] = hy[mm] + (ez[mm + 1] - ez[mm]) / imp0;

    // update electric field
    for (auto mm{ 1 }; mm < size; ++mm)
      ez[mm] = ez[mm] + (hy[mm] - hy[mm - 1]) * imp0;

    // source
    ez[0] = exp(-(time - 30.0) * (time - 30.0) / 100.0);

    std::copy(ez.begin(), ez.end(), history_ez.begin() + time * size);
  }

  vtkNew<vtkImageData> image;

  image->SetDimensions(size, size, 1);
  // image->SetSpacing(1.0, 1.0, 1.0);
  // image->SetOrigin(0.0, 0.0, 0.0);
  image->AllocateScalars(VTK_FLOAT, 3);

  mkdir("out", 0755);
  mkdir("out/animation", 0755);

  for (auto time{ 0 }; time < max_time; ++time) {
    const int* dims = image->GetDimensions();

    for (auto y = 0; y < dims[1]; ++y) {
      for (auto x = 0; x < dims[0]; ++x) {
        float* pixel = static_cast<float*>(image->GetScalarPointer(x, y, 0));
        pixel[0] = 0.0;
        pixel[1] = 0.0;
        pixel[2] = 0.0;
      }
    }

    for (auto x{ 0 }; x < dims[0]; ++x) {
      float* pixel = static_cast<float*>(image->GetScalarPointer(x, 0, 0));
      pixel[0] = 0.0;
      pixel[1] = history_ez[size * time + x];
      pixel[2] = 0.0;
    }

    image->GetPointData()->SetActiveVectors("ImageScalars");

    image->Modified();

    vtkNew<vtkXMLImageDataWriter> writer;
    const std::string file_name = "out/animation/fdtd" + std::to_string(time) + ".vti";
    writer->SetFileName(file_name.c_str());
    writer->SetInputData(image);
    writer->Write();
  }

  std::ofstream output;
  output.open("out/fdtd.pvd");

  output << "<?xml version=\"1.0\"?>" << std::endl;
  output << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  output << "  <Collection>" << std::endl;

  for (auto t{ 0 }; t < max_time; ++t) {
    output << "    <DataSet timestep=\"" << std::to_string(t) << "\" group=\"\" part=\"0\" file=\"animation/fdtd" << std::to_string(t) << ".vti\" name=\"fdtd\"/>" << std::endl;
  }

  output << "  </Collection>" << std::endl;
  output << "</VTKFile>" << std::endl;

  output.close();

  return 0;
}


