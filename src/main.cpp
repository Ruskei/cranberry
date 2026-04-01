#include <chrono>
#include <iostream>
#include <ostream>
#include <sys/stat.h>
#include <sys/types.h>

#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>

#include "config.hpp"
#include "grid.hpp"
#include "particles.hpp"
#include "print_result.hpp"
#include "field_view_2d_result.hpp"

void run_3d() {
  std::cout << "Running 3D FDTD sim!!" << std::endl;

  const int size = Config::size;

  std::vector<SimulationResult<size>> results;
  results.emplace_back(PrintResult{10});
  results.emplace_back(FieldView2D<size>("e-field", WhichField::E, Axis::Z, 64, 10));

  Sim<size> grid(
    std::vector<Particle>{
      Particle{64, 16, 64, 0, 0.5, 0, 1, 1},
      // Particle{64, 70, 64, 0, -0.5, 0, -1, 1},
    },
    std::move(results)
  );

  grid.setup_coefficients();

  auto start = std::chrono::steady_clock::now();

  while (grid.time < Config::max_time)
    grid.step();
  
  grid.finish();

  auto finished = std::chrono::steady_clock::now();

  std::cout << "Sim finished!" << std::endl;
  std::cout << "took: " << std::chrono::duration_cast<std::chrono::milliseconds>(finished - start) << "ms" << std::endl;
}

int main() {
  run_3d();

  return 0;
}
