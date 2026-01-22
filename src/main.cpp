#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
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
constexpr int max_time = 300;

namespace {
void run_3d();
}

int main() {
  std::cout << "Running simulation\n";

  std::vector<double> history_ez(size * max_time);

  std::array<double, size> ez{};
  std::array<double, size> hy{};

  for (auto time{0}; time < max_time; ++time) {
    hy[size - 1] = hy[size - 2];
    // update magnetic field
    for (auto mm{0}; mm < size - 1; ++mm)
      hy[mm] = hy[mm] + (ez[mm + 1] - ez[mm]) / imp0;

    hy[49] -= exp(-(time - 30.0) * (time - 30.0) / 100.0) / imp0;
    ez[0] = ez[1];

    // update electric field
    for (auto mm{1}; mm < size; ++mm)
      ez[mm] = ez[mm] + (hy[mm] - hy[mm - 1]) * imp0;

    // source
    ez[50] += exp(-(time + 0.5 - (-0.5) - 30.0) * (time + 0.5 - (0.5) - 30.0) /
                  100.0);

    std::copy(ez.begin(), ez.end(), history_ez.begin() + time * size);
  }

  vtkNew<vtkImageData> image;

  image->SetDimensions(size, size, 1);
  // image->SetSpacing(1.0, 1.0, 1.0);
  // image->SetOrigin(0.0, 0.0, 0.0);
  image->AllocateScalars(VTK_FLOAT, 3);

  mkdir("out", 0755);
  mkdir("out/animation", 0755);

  for (auto time{0}; time < max_time; ++time) {
    const int *dims = image->GetDimensions();

    for (auto y = 0; y < dims[1]; ++y) {
      for (auto x = 0; x < dims[0]; ++x) {
        float *pixel = static_cast<float *>(image->GetScalarPointer(x, y, 0));
        pixel[0] = 0.0;
        pixel[1] = 0.0;
        pixel[2] = 0.0;
      }
    }

    for (auto x{0}; x < dims[0]; ++x) {
      float *pixel = static_cast<float *>(image->GetScalarPointer(x, 0, 0));
      pixel[0] = history_ez[size * time + x];
      pixel[1] = history_ez[size * time + x];
      pixel[2] = 0.0;
    }

    image->GetPointData()->SetActiveVectors("ImageScalars");

    image->Modified();

    vtkNew<vtkXMLImageDataWriter> writer;
    const std::string file_name =
        "out/animation/fdtd" + std::to_string(time) + ".vti";
    writer->SetFileName(file_name.c_str());
    writer->SetInputData(image);
    writer->Write();
  }

  std::ofstream output;
  output.open("out/fdtd.pvd");

  output << "<?xml version=\"1.0\"?>" << std::endl;
  output << "<VTKFile type=\"Collection\" version=\"0.1\" "
            "byte_order=\"LittleEndian\">"
         << std::endl;
  output << "  <Collection>" << std::endl;

  for (auto t{0}; t < max_time; ++t) {
    output << "    <DataSet timestep=\"" << std::to_string(t)
           << "\" group=\"\" part=\"0\" file=\"animation/fdtd"
           << std::to_string(t) << ".vti\" name=\"fdtd\"/>" << std::endl;
  }

  output << "  </Collection>" << std::endl;
  output << "</VTKFile>" << std::endl;

  output.close();

  run_3d();

  return 0;
}

constexpr int size_3d = 40;
const double courants = 1.0 / sqrt(3.0);

constexpr double f0 = 0.08;
constexpr double t0 = 30.0;

int idx_hx(int x, int y, int z) {
  return (x * (size_3d - 1) + y) * (size_3d - 1) + z;
}
int idx_hy(int x, int y, int z) {
  return (x * size_3d + y) * (size_3d - 1) + z;
}
int idx_hz(int x, int y, int z) {
  return (x * (size_3d - 1) + y) * size_3d + z;
}

int idx_ex(int x, int y, int z) { return (x * size_3d + y) * size_3d + z; }
int idx_ey(int x, int y, int z) {
  return (x * (size_3d - 1) + y) * size_3d + z;
}
int idx_ez(int x, int y, int z) {
  return (x * size_3d + y) * (size_3d - 1) + z;
}

// abc
int idx_eyx(int a, int b) { return a * size_3d + b; }
int idx_ezx(int a, int b) { return a * (size_3d - 1) + b; }

int idx_exy(int a, int b) { return a * size_3d + b; }
int idx_ezy(int a, int b) { return a * (size_3d - 1) + b; }

int idx_exz(int a, int b) { return a * size_3d + b; }
int idx_eyz(int a, int b) { return a * (size_3d - 1) + b; }

const double abcco = (courants - 1.0) / (courants + 1.0);

namespace {
void run_3d() {
  std::cout << "Running 3D FDTD sim!!" << std::endl;

  constexpr int e_size = (size_3d - 1) * size_3d * size_3d;
  constexpr int h_size = (size_3d - 1) * (size_3d - 1) * size_3d;

  std::vector<double> history_ex(max_time * e_size);
  std::vector<double> history_ey(max_time * e_size);
  std::vector<double> history_ez(max_time * e_size);

  std::vector<double> history_hx(max_time * h_size);
  std::vector<double> history_hy(max_time * h_size);
  std::vector<double> history_hz(max_time * h_size);

  std::vector<double> ex(e_size);
  std::vector<double> cexe(e_size);
  std::vector<double> cexh(e_size);

  std::vector<double> ey(e_size);
  std::vector<double> ceye(e_size);
  std::vector<double> ceyh(e_size);

  std::vector<double> ez(e_size);
  std::vector<double> ceze(e_size);
  std::vector<double> cezh(e_size);

  std::vector<double> hx(h_size);
  std::vector<double> chxh(h_size);
  std::vector<double> chxe(h_size);

  std::vector<double> hy(h_size);
  std::vector<double> chyh(h_size);
  std::vector<double> chye(h_size);

  std::vector<double> hz(h_size);
  std::vector<double> chzh(h_size);
  std::vector<double> chze(h_size);

  // set up coefficients
  for (auto x{0}; x < size_3d - 1; ++x)
    for (auto y{0}; y < size_3d; ++y)
      for (auto z{0}; z < size_3d; ++z) {
        const int i = idx_ex(x, y, z);
        cexe[i] = 1.0;
        cexh[i] = courants * imp0;
      }

  for (auto x{0}; x < size_3d; ++x)
    for (auto y{0}; y < size_3d - 1; ++y)
      for (auto z{0}; z < size_3d; ++z) {
        const int i = idx_ey(x, y, z);
        ceye[i] = 1.0;
        ceyh[i] = courants * imp0;
      }

  for (auto x{0}; x < size_3d; ++x)
    for (auto y{0}; y < size_3d; ++y)
      for (auto z{0}; z < size_3d - 1; ++z) {
        const int i = idx_ez(x, y, z);
        ceze[i] = 1.0;
        cezh[i] = courants * imp0;
      }

  for (auto x{0}; x < size_3d; ++x)
    for (auto y{0}; y < size_3d - 1; ++y)
      for (auto z{0}; z < size_3d - 1; ++z) {
        const int i = idx_hx(x, y, z);
        chxh[i] = 1.0;
        chxe[i] = courants / imp0;
      }

  for (auto x{0}; x < size_3d - 1; ++x)
    for (auto y{0}; y < size_3d; ++y)
      for (auto z{0}; z < size_3d - 1; ++z) {
        const int i = idx_hy(x, y, z);
        chyh[i] = 1.0;
        chye[i] = courants / imp0;
      }

  for (auto x{0}; x < size_3d - 1; ++x)
    for (auto y{0}; y < size_3d - 1; ++y)
      for (auto z{0}; z < size_3d; ++z) {
        const int i = idx_hz(x, y, z);
        chzh[i] = 1.0;
        chze[i] = courants / imp0;
      }

  const int size_abc = (size_3d - 1) * size_3d;

  std::vector<double> eyx0(size_abc);
  std::vector<double> ezx0(size_abc);
  std::vector<double> eyx1(size_abc);
  std::vector<double> ezx1(size_abc);

  std::vector<double> exy0(size_abc);
  std::vector<double> ezy0(size_abc);
  std::vector<double> exy1(size_abc);
  std::vector<double> ezy1(size_abc);

  std::vector<double> exz0(size_abc);
  std::vector<double> eyz0(size_abc);
  std::vector<double> exz1(size_abc);
  std::vector<double> eyz1(size_abc);

  for (auto time{0}; time < max_time; ++time) {
    const double t = time * courants;
    // update magnetic fields
    for (auto x{0}; x < size_3d; ++x)
      for (auto y{0}; y < size_3d - 1; ++y)
        for (auto z{0}; z < size_3d - 1; ++z) {
          int i = idx_hx(x, y, z);
          hx[i] = chxh[i] * hx[i] +
                  chxe[i] * ((ey[idx_ey(x, y, z + 1)] - ey[idx_ey(x, y, z)]) -
                             (ez[idx_ez(x, y + 1, z)] - ez[idx_ez(x, y, z)]));
        }

    for (auto x{0}; x < size_3d - 1; ++x)
      for (auto y{0}; y < size_3d; ++y)
        for (auto z{0}; z < size_3d - 1; ++z) {
          int i = idx_hy(x, y, z);
          hy[i] = chyh[i] * hy[i] +
                  chye[i] * ((ez[idx_ez(x + 1, y, z)] - ez[idx_ez(x, y, z)]) -
                             (ex[idx_ex(x, y, z + 1)] - ex[idx_ex(x, y, z)]));
        }

    for (auto x{0}; x < size_3d - 1; ++x)
      for (auto y{0}; y < size_3d - 1; ++y)
        for (auto z{0}; z < size_3d; ++z) {
          int i = idx_hz(x, y, z);
          hz[i] = chzh[i] * hz[i] +
                  chze[i] * ((ex[idx_ex(x, y + 1, z)] - ex[idx_ex(x, y, z)]) -
                             (ey[idx_ey(x + 1, y, z)] - ey[idx_ey(x, y, z)]));
        }

    // update electric fields
    for (auto x{0}; x < size_3d - 1; ++x)
      for (auto y{1}; y < size_3d - 1; ++y)
        for (auto z{1}; z < size_3d - 1; ++z) {
          int i = idx_ex(x, y, z);
          ex[i] = cexe[i] * ex[i] +
                  cexh[i] * ((hz[idx_hz(x, y, z)] - hz[idx_hz(x, y - 1, z)]) -
                             (hy[idx_hy(x, y, z)] - hy[idx_hy(x, y, z - 1)]));
        }

    for (auto x{1}; x < size_3d - 1; ++x)
      for (auto y{0}; y < size_3d - 1; ++y)
        for (auto z{1}; z < size_3d - 1; ++z) {
          int i = idx_ey(x, y, z);
          ey[i] = ceye[i] * ey[i] +
                  ceyh[i] * ((hx[idx_hx(x, y, z)] - hx[idx_hx(x, y, z - 1)]) -
                             (hz[idx_hz(x, y, z)] - hz[idx_hz(x - 1, y, z)]));
        }

    for (auto x{1}; x < size_3d - 1; ++x)
      for (auto y{1}; y < size_3d - 1; ++y)
        for (auto z{0}; z < size_3d - 1; ++z) {
          int i = idx_ez(x, y, z);
          ez[i] = ceze[i] * ez[i] +
                  cezh[i] * ((hy[idx_hy(x, y, z)] - hy[idx_hy(x - 1, y, z)]) -
                             (hx[idx_hx(x, y, z)] - hx[idx_hx(x, y - 1, z)]));
        }

    const int insertion_idx =
        idx_ex((size_3d - 1) / 2, (size_3d - 1) / 2, (size_3d - 1) / 2);

    double tau = t - t0;
    double a = M_PI * M_PI * f0 * f0 * tau * tau;
    ex[insertion_idx] += (1.0 - 2.0 * a) * std::exp(-a);

    // pec
    // for (auto y{0}; y < size_3d - 1; ++y)
    //   for (auto z{0}; z < size_3d; ++z) {
    //     ey[idx_ey(0, y, z)] = 0.0;
    //     ey[idx_ey(size_3d - 1, y, z)] = 0.0;
    //   }
    //
    // for (auto y{0}; y < size_3d; ++y)
    //   for (auto z{0}; z < size_3d - 1; ++z) {
    //     ez[idx_ez(0, y, z)] = 0.0;
    //     ez[idx_ez(size_3d - 1, y, z)] = 0.0;
    //   }
    //
    // for (auto x{0}; x < size_3d - 1; ++x)
    //   for (auto z{0}; z < size_3d; ++z) {
    //     ex[idx_ex(x, 0, z)] = 0.0;
    //     ex[idx_ex(x, size_3d - 1, z)] = 0.0;
    //   }
    //
    // for (auto x{0}; x < size_3d; ++x)
    //   for (auto z{0}; z < size_3d - 1; ++z) {
    //     ez[idx_ez(x, 0, z)] = 0.0;
    //     ez[idx_ez(x, size_3d - 1, z)] = 0.0;
    //   }
    //
    // for (auto x{0}; x < size_3d - 1; ++x)
    //   for (auto y{0}; y < size_3d; ++y) {
    //     ex[idx_ex(x, y, 0)] = 0.0;
    //     ex[idx_ex(x, y, size_3d - 1)] = 0.0;
    //   }
    //
    // for (auto x{0}; x < size_3d; ++x)
    //   for (auto y{0}; y < size_3d - 1; ++y) {
    //     ey[idx_ey(x, y, 0)] = 0.0;
    //     ey[idx_ey(x, y, size_3d - 1)] = 0.0;
    //   }

    // abc
    /* ABC at "x0" */
    {
      const int x{0};
      for (int y{0}; y < size_3d - 1; ++y)
        for (int z{0}; z < size_3d; ++z) {
          ey[idx_ey(x, y, z)] =
              eyx0[idx_eyx(y, z)] +
              abcco * (ey[idx_ey(x + 1, y, z)] - ey[idx_ey(x, y, z)]);
          eyx0[idx_eyx(y, z)] = ey[idx_ey(x + 1, y, z)];
        }

      for (int y{0}; y < size_3d; ++y)
        for (int z{0}; z < size_3d - 1; ++z) {
          ez[idx_ez(x, y, z)] =
              ezx0[idx_ezx(y, z)] +
              abcco * (ez[idx_ez(x + 1, y, z)] - ez[idx_ez(x, y, z)]);
          ezx0[idx_ezx(y, z)] = ez[idx_ez(x + 1, y, z)];
        }
    }

    /* ABC at "x1" */
    {
      const int x = size_3d - 1;
      for (int y{0}; y < size_3d - 1; ++y)
        for (int z{0}; z < size_3d; ++z) {
          ey[idx_ey(x, y, z)] =
              eyx1[idx_eyx(y, z)] +
              abcco * (ey[idx_ey(x - 1, y, z)] - ey[idx_ey(x, y, z)]);
          eyx1[idx_eyx(y, z)] = ey[idx_ey(x - 1, y, z)];
        }

      for (int y{0}; y < size_3d; ++y)
        for (int z{0}; z < size_3d - 1; ++z) {
          ez[idx_ez(x, y, z)] =
              ezx1[idx_ezx(y, z)] +
              abcco * (ez[idx_ez(x - 1, y, z)] - ez[idx_ez(x, y, z)]);
          ezx1[idx_ezx(y, z)] = ez[idx_ez(x - 1, y, z)];
        }
    }

    /* ABC at "y0" */
    {
      const int y{0};
      for (int x{0}; x < size_3d - 1; ++x)
        for (int z{0}; z < size_3d; ++z) {
          ex[idx_ex(x, y, z)] =
              exy0[idx_exy(x, z)] +
              abcco * (ex[idx_ex(x, y + 1, z)] - ex[idx_ex(x, y, z)]);
          exy0[idx_exy(x, z)] = ex[idx_ex(x, y + 1, z)];
        }

      for (int x{0}; x < size_3d; ++x)
        for (int z{0}; z < size_3d - 1; ++z) {
          ez[idx_ez(x, y, z)] =
              ezy0[idx_ezy(x, z)] +
              abcco * (ez[idx_ez(x, y + 1, z)] - ez[idx_ez(x, y, z)]);
          ezy0[idx_ezy(x, z)] = ez[idx_ez(x, y + 1, z)];
        }
    }

    /* ABC at "y1" */
    {
      const int y = size_3d - 1;
      for (int x{0}; x < size_3d - 1; ++x)
        for (int z{0}; z < size_3d; ++z) {
          ex[idx_ex(x, y, z)] =
              exy1[idx_exy(x, z)] +
              abcco * (ex[idx_ex(x, y - 1, z)] - ex[idx_ex(x, y, z)]);
          exy1[idx_exy(x, z)] = ex[idx_ex(x, y - 1, z)];
        }

      for (int x{0}; x < size_3d; ++x)
        for (int z{0}; z < size_3d - 1; ++z) {
          ez[idx_ez(x, y, z)] =
              ezy1[idx_ezy(x, z)] +
              abcco * (ez[idx_ez(x, y - 1, z)] - ez[idx_ez(x, y, z)]);
          ezy1[idx_ezy(x, z)] = ez[idx_ez(x, y - 1, z)];
        }
    }

    /* ABC at "z0" (bottom) */
    {
      const int z{0};
      for (int x{0}; x < size_3d - 1; ++x)
        for (int y{0}; y < size_3d; ++y) {
          ex[idx_ex(x, y, z)] =
              exz0[idx_exz(x, y)] +
              abcco * (ex[idx_ex(x, y, z + 1)] - ex[idx_ex(x, y, z)]);
          exz0[idx_exz(x, y)] = ex[idx_ex(x, y, z + 1)];
        }

      for (int x{0}; x < size_3d; ++x)
        for (int y{0}; y < size_3d - 1; ++y) {
          ey[idx_ey(x, y, z)] =
              eyz0[idx_eyz(x, y)] +
              abcco * (ey[idx_ey(x, y, z + 1)] - ey[idx_ey(x, y, z)]);
          eyz0[idx_eyz(x, y)] = ey[idx_ey(x, y, z + 1)];
        }
    }

    /* ABC at "z1" (top) */
    {
      const int z = size_3d - 1;
      for (int x{0}; x < size_3d - 1; ++x)
        for (int y{0}; y < size_3d; ++y) {
          ex[idx_ex(x, y, z)] =
              exz1[idx_exz(x, y)] +
              abcco * (ex[idx_ex(x, y, z - 1)] - ex[idx_ex(x, y, z)]);
          exz1[idx_exz(x, y)] = ex[idx_ex(x, y, z - 1)];
        }

      for (int x{0}; x < size_3d; ++x)
        for (int y{0}; y < size_3d - 1; ++y) {
          ey[idx_ey(x, y, z)] =
              eyz1[idx_eyz(x, y)] +
              abcco * (ey[idx_ey(x, y, z - 1)] - ey[idx_ey(x, y, z)]);
          eyz1[idx_eyz(x, y)] = ey[idx_ey(x, y, z - 1)];
        }
    }

    std::copy(ex.begin(), ex.end(), history_ex.begin() + time * e_size);
    std::copy(ey.begin(), ey.end(), history_ey.begin() + time * e_size);
    std::copy(ez.begin(), ez.end(), history_ez.begin() + time * e_size);

    std::copy(hx.begin(), hx.end(), history_hx.begin() + time * h_size);
    std::copy(hy.begin(), hy.end(), history_hy.begin() + time * h_size);
    std::copy(hz.begin(), hz.end(), history_hz.begin() + time * h_size);
  }

  std::cout << "Finished 3D simulation!" << std::endl;

  for (auto time{0}; time < max_time; ++time) {
    vtkNew<vtkImageData> image;

    image->SetDimensions(size_3d, size_3d, size_3d);
    image->AllocateScalars(VTK_FLOAT, 3);

    mkdir("out", 0755);
    mkdir("out/3d_animation", 0755);

    const int *dims = image->GetDimensions();

    for (auto z{0}; z < dims[2]; ++z)
      for (auto y{0}; y < dims[1]; ++y)
        for (auto x{0}; x < dims[0]; ++x) {
          float *pixel = static_cast<float *>(image->GetScalarPointer(x, y, z));
          pixel[0] = 0;
          pixel[1] = 0;
          pixel[2] = 0;
        }

    for (auto z{1}; z < dims[2] - 1; ++z)
      for (auto y{1}; y < dims[1] - 1; ++y)
        for (auto x{1}; x < dims[0] - 1; ++x) {
          float *pixel = static_cast<float *>(image->GetScalarPointer(x, y, z));
          pixel[0] = history_ex[e_size * time + idx_ex(x, y, z)];
          pixel[1] = history_ey[e_size * time + idx_ey(x, y, z)];
          pixel[2] = history_ez[e_size * time + idx_ez(x, y, z)];
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

  std::ofstream output;
  output.open("out/3d_animation.pvd");

  output << "<?xml version=\"1.0\"?>" << std::endl;
  output << "<VTKFile type=\"Collection\" version=\"0.1\" "
            "byte_order=\"LittleEndian\">"
         << std::endl;
  output << "  <Collection>" << std::endl;

  for (auto t{0}; t < max_time; ++t) {
    output << "    <DataSet timestep=\"" << std::to_string(t)
           << "\" group=\"\" part=\"0\" file=\"3d_animation/efield"
           << std::to_string(t) << ".vti\" name=\"efield\"/>" << std::endl;
  }

  output << "  </Collection>" << std::endl;
  output << "</VTKFile>" << std::endl;

  output.close();
}
} // namespace
