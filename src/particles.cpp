#include "particles.hpp"

#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>

void ParticleWriter::write_particles(const std::vector<Particle> &particles,
                                     double time) {
  vtkNew<vtkPoints> points;
  vtkNew<vtkFloatArray> velocity;
  velocity->SetName("velocity");
  velocity->SetNumberOfComponents(3);
  for (const auto &p : particles) {
    points->InsertNextPoint(p.px, p.py, p.pz);
    velocity->InsertNextTuple3(p.vx, p.vy, p.vz);
  }

  vtkNew<vtkPolyData> poly;
  poly->SetPoints(points);
  poly->GetPointData()->AddArray(velocity);
  poly->GetPointData()->SetActiveVectors("velocity");

  vtkNew<vtkXMLPolyDataWriter> writer;
  const std::string filename =
      "out/" + name + "/particles" + std::to_string(time) + ".vtp";
  writer->SetFileName(filename.c_str());
  writer->SetInputData(poly);
  writer->Write();
}

void ParticleWriter::write_all(double max_time, double dt) {
  std::ofstream output;
  output.open(("out/" + name + ".pvd").c_str());

  output << "<?xml version=\"1.0\"?>" << std::endl;
  output << "<VTKFile type=\"Collection\" version=\"0.1\" "
            "byte_order=\"LittleEndian\">"
         << std::endl;
  output << "  <Collection>" << std::endl;

  for (auto t{0.0}; t < max_time; t += dt) {
    output << "    <DataSet timestep=\"" << std::to_string(t)
           << "\" group=\"\" part=\"0\" file=\"" + name + "/particles"
           << std::to_string(t) << ".vtp\" name=\"particles\"/>" << std::endl;
  }

  output << "  </Collection>" << std::endl;
  output << "</VTKFile>" << std::endl;

  output.close();
}

double form_factor(double px, int nx) {
  const int i = static_cast<int>(std::floor(px));

  if (nx < i || nx > i + 1)
    return 0.0;

  const double fx = px - i;

  return (i == nx) ? 1.0 - fx : fx;
}

double form_factor_diff_helper(double s_o_x, double s_o_y, double s_o_z,
                               double s_n_x, double s_n_y, double s_n_z) {
  return (s_n_x * s_n_y * s_n_z - s_o_x * s_n_y * s_n_z +
          s_n_x * s_o_y * s_o_z - s_o_x * s_o_y * s_o_z) /
             3.0 +
         (s_n_x * s_o_y * s_n_z - s_o_x * s_o_y * s_n_z +
          s_n_x * s_n_y * s_o_z - s_o_x * s_n_y * s_o_z) /
             6.0;
}

