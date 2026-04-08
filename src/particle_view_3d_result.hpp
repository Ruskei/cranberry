#pragma once

#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <utility>
#include <vector>

#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>

#include "config.hpp"
#include "particles.hpp"

enum class ParticleVelocityOutput { Coordinate, Proper };

struct ParticleView3D {
  std::string name;
  int interval{1};
  int num_saved{};
  ParticleVelocityOutput velocity_output{ParticleVelocityOutput::Coordinate};

  ParticleView3D(std::string name, int interval,
                 ParticleVelocityOutput velocity_output =
                     ParticleVelocityOutput::Coordinate)
      : name{std::move(name)}, interval{interval},
        velocity_output{velocity_output} {
    mkdir("out", 0755);
    mkdir(("out/" + this->name).c_str(), 0755);
    mkdir(("out/" + this->name + "/particles").c_str(), 0755);
  }

  std::array<double, 3> output_velocity(const Particle &particle) const {
    if (velocity_output == ParticleVelocityOutput::Coordinate)
      return {particle.vx, particle.vy, particle.vz};

    return {particle.ux, particle.uy, particle.uz};
  }

  void step(double time, const std::vector<Particle> &particles) {
    if (static_cast<int>(std::round(time / Config::dt)) % interval != 0)
      return;

    vtkNew<vtkPoints> points;
    vtkNew<vtkFloatArray> velocity;
    velocity->SetName("velocity");
    velocity->SetNumberOfComponents(3);

    for (const auto &particle : particles) {
      points->InsertNextPoint(particle.px, particle.py, particle.pz);
      const auto v = output_velocity(particle);
      velocity->InsertNextTuple3(v[0], v[1], v[2]);
    }

    vtkNew<vtkPolyData> poly;
    poly->SetPoints(points);
    poly->GetPointData()->AddArray(velocity);
    poly->GetPointData()->SetActiveVectors("velocity");

    vtkNew<vtkXMLPolyDataWriter> writer;
    const std::string file_name =
        "out/" + name + "/particles/" + std::to_string(num_saved++) + ".vtp";
    writer->SetFileName(file_name.c_str());
    writer->SetInputData(poly);
    writer->Write();
  }

  void finish() {
    std::ofstream output(("out/" + name + "/recording.pvd").c_str());

    output << "<?xml version=\"1.0\"?>" << std::endl;
    output << "<VTKFile type=\"Collection\" version=\"0.1\" "
              "byte_order=\"LittleEndian\">"
           << std::endl;
    output << "  <Collection>" << std::endl;

    for (auto i{0}; i < num_saved; ++i) {
      output << "    <DataSet timestep=\"" << std::to_string(i)
             << "\" group=\"\" part=\"0\" file=\"particles/"
             << std::to_string(i) << ".vtp\" name=\"particle_view_3d\"/>"
             << std::endl;
    }

    output << "  </Collection>" << std::endl;
    output << "</VTKFile>" << std::endl;
  }
};
