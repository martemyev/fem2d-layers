#include "result.h"
#include "auxiliary_functions.h"
#include "dof_handler.h"
#include "fine_mesh.h"
#include "point.h"
#include <boost/filesystem.hpp>
#include <fstream>


Result::Result(const DoFHandler *dof_handler)
  : _dof_handler(dof_handler)
{ }



Result::~Result()
{ }



void Result::write_vtu(const std::string &filename,
                       const Vec &solution,
                       const Vec &exact_solution) const
{
  using namespace boost::filesystem;
  expect(extension(filename) == ".vtu", "The extension of the file ('" + extension(filename) + "') is not suitable for this function.");

  expect(_dof_handler->n_dofs() == _dof_handler->fmesh()->n_vertices(),
         "This function should be rewritten for the case higher order basis functions"
         "(when the number of degrees of freedom is not equal to the number of the mesh vertices");

  std::ofstream out(filename.c_str()); // open the file for writing
  require(out, "File " + filename + " cannot be opened");

  // extract the data from PETSc vectors
  std::vector<int> idx(_dof_handler->n_dofs());
  std::iota(idx.begin(), idx.end(), 0); // idx = { 0, 1, 2, 3, .... }
  std::vector<double> solution_values(_dof_handler->n_dofs());
  VecGetValues(solution, _dof_handler->n_dofs(), &idx[0], &solution_values[0]);

  std::vector<double> exact_solution_values;
  if (exact_solution)
  {
    exact_solution_values.resize(_dof_handler->n_dofs());
    VecGetValues(exact_solution, _dof_handler->n_dofs(), &idx[0], &exact_solution_values[0]);
  }

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  out << "  <UnstructuredGrid>\n";
  out << "    <Piece NumberOfPoints=\"" << _dof_handler->fmesh()->n_vertices() << "\" NumberOfCells=\"" << _dof_handler->fmesh()->n_triangles() << "\">\n";
  out << "      <PointData>\n";
  out << "        <DataArray type=\"Float64\" Name=\"solution\" format=\"ascii\">\n";
  for (int i = 0; i < _dof_handler->fmesh()->n_vertices(); ++i)
    out << solution_values[i] << "\n";
  out << "        </DataArray>\n";
  if (exact_solution)
  {
    out << "        <DataArray type=\"Float64\" Name=\"exact_solution\" format=\"ascii\">\n";
    for (int i = 0; i < _dof_handler->fmesh()->n_vertices(); ++i)
      out << exact_solution_values[i] << "\n";
    out << "        </DataArray>\n";
  }
  out << "      </PointData>\n";
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"" << Point::n_coord << "\" format=\"ascii\">\n";
  const std::vector<Point> vertices = _dof_handler->fmesh()->vertices();
  for (int i = 0; i < _dof_handler->fmesh()->n_vertices(); ++i)
  {
    for (int j = 0; j < Point::n_coord; ++j)
      out << vertices[i].coord(j) << " ";
    out << "\n";
  }
  out << "        </DataArray>\n";
  out << "      </Points>\n";
  out << "      <Cells>\n";
  out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for (int tr = 0; tr < _dof_handler->fmesh()->n_triangles(); ++tr)
  {
    for (int vert = 0; vert < Triangle::n_vertices; ++vert)
      out << _dof_handler->fmesh()->triangle(tr).vertex(vert) << " ";
    out << "\n";
  }
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for (int tr = 0; tr < _dof_handler->fmesh()->n_triangles(); ++tr)
    out << (tr + 1) * Triangle::n_vertices << " ";
  out << "\n";
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
  for (int tr = 0; tr < _dof_handler->fmesh()->n_triangles(); ++tr)
    out << Triangle::vtk_el_type << " ";
  out << "\n";
  out << "        </DataArray>\n";
  out << "      </Cells>\n";
  out << "    </Piece>\n";
  out << "  </UnstructuredGrid>\n";
  out << "</VTKFile>\n";

  out.close();
}
