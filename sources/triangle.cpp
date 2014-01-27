#include "triangle.h"
#include "auxiliary_functions.h"
#include <cmath>


Triangle::Triangle()
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type),
    _detD(0.)
{
  for (int i = 0; i < n_vertices; ++i)
    _A[i] = _B[i] = _C[i] = 0.;
}



Triangle::~Triangle()
{
  _ghost_cells.clear();
  _dofs.clear();
}



Triangle::Triangle(const std::vector<unsigned int> &ver,
                   const std::vector<Point> &mesh_vertices,
                   const unsigned int mat_id,
                   const unsigned int part_id,
                   const std::vector<unsigned int> &g_cells)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  expect(_vertices.size() == ver.size(),
         "The size of list of vertices (" + d2s(ver.size()) +
         ") is not equal to really needed number of vertices (" + d2s(n_vertices) + ")");
  _vertices = ver;
  _material_id = mat_id;
  _partition_id = part_id;
  _ghost_cells = g_cells;

  if (!mesh_vertices.empty())
    calculate_barycentric(mesh_vertices);
}



Triangle::Triangle(const unsigned int v1,
                   const unsigned int v2,
                   const unsigned int v3,
                   const std::vector<Point> &mesh_vertices,
                   const unsigned int mat_id,
                   const unsigned int part_id,
                   const std::vector<unsigned int> &g_cells)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  _vertices[0] = v1;
  _vertices[1] = v2;
  _vertices[2] = v3;
  _material_id = mat_id;
  _partition_id = part_id;
  _ghost_cells = g_cells;

  if (!mesh_vertices.empty())
    calculate_barycentric(mesh_vertices);
}



Triangle::Triangle(const Triangle &tri)
  : MeshElement(tri),
    _ghost_cells(tri._ghost_cells),
    _dofs(tri._dofs),
    _detD(tri._detD)
{
  for (int i = 0; i < n_vertices; ++i)
  {
    _A[i] = tri._A[i];
    _B[i] = tri._B[i];
    _C[i] = tri._C[i];
  }
}



Triangle& Triangle::operator =(const Triangle &tri)
{
  MeshElement::operator =(tri);
  _ghost_cells = tri._ghost_cells;
  _dofs = tri._dofs;
  _detD = tri._detD;
  for (int i = 0; i < n_vertices; ++i)
  {
    _A[i] = tri._A[i];
    _B[i] = tri._B[i];
    _C[i] = tri._C[i];
  }
}



unsigned int Triangle::n_ghost_cells() const
{
  return _ghost_cells.size();
}


unsigned int Triangle::ghost_cell(unsigned int num) const
{
  expect(num >= 0 && num < _ghost_cells.size(), "Incorrect input parameter");
  return _ghost_cells[num];
}



unsigned int Triangle::n_dofs() const
{
  return _dofs.size();
}



unsigned int Triangle::dof(unsigned int num) const
{
  expect(num >= 0 && num < _dofs.size(), "Incorrect input parameter");
  return _dofs[num];
}



void Triangle::n_dofs(unsigned int number)
{
  _dofs.resize(number);
}



void Triangle::dof(unsigned int num, unsigned int value)
{
  expect(num >= 0 && num < _dofs.size(), "Incorrect input parameter");
  _dofs[num] = value;
}



void Triangle::local_mass_matrix(double **loc_mat) const
{
  expect(fabs(_detD) > 1e-15,
         "An attempt to calculate a local matrix for a singular triangle (detD = " + d2s(_detD, true) + ")");

  switch (_dofs.size())
  {
  case n_dofs_first:
  {
    const double mat[][n_dofs_first] = { { 2., 1., 1. },
                                         { 1., 2., 1. },
                                         { 1., 1., 2. }
                                       };
    const double rho = 1.; // averaged coefficient
    for (int i = 0; i < n_dofs_first; ++i)
      for (int j = 0; j < n_dofs_first; ++j)
        loc_mat[i][j] = rho * fabs(_detD) * mat[i][j] / 24.;
    return;
  }
  default:
    require(false, "Unknown fe order for generating local matrix");
  }
}



void Triangle::local_stiffness_matrix(double **loc_mat, double coef_a) const
{
  expect(fabs(_detD) > 1e-15,
         "An attempt to calculate a local matrix for a singular triangle (detD = " + d2s(_detD, true) + ")");

  switch (_dofs.size())
  {
  case n_dofs_first:
  {
    for (int i = 0; i < n_dofs_first; ++i)
      for (int j = 0; j < n_dofs_first; ++j)
        loc_mat[i][j] = coef_a * fabs(_detD) * (_A[i]*_A[j] + _B[i]*_B[j]) / 2.;
    return;
  }
  default:
    require(false, "Unknown fe order for generating local matrix");
  }
}



void Triangle::calculate_barycentric(const std::vector<Point> &mesh_vertices)
{
  double x[n_vertices], y[n_vertices]; // coordinates of the triangle's vertices
  for (int i = 0; i < n_vertices; ++i)
  {
    const unsigned int vertex = _vertices[i];
    x[i] = mesh_vertices[vertex].coord(0); // x-coordinate
    y[i] = mesh_vertices[vertex].coord(1); // y-coordinate
  }
  _detD = (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]);
  _A[0] = (y[1] - y[2]) / _detD;
  _A[1] = (y[2] - y[0]) / _detD;
  _A[2] = (y[0] - y[1]) / _detD;
  _B[0] = (x[2] - x[1]) / _detD;
  _B[1] = (x[0] - x[2]) / _detD;
  _B[2] = (x[1] - x[0]) / _detD;
  _C[0] = (x[1]*y[2] - x[2]*y[1]) / _detD;
  _C[1] = (x[2]*y[0] - x[0]*y[2]) / _detD;
  _C[2] = (x[0]*y[1] - x[1]*y[0]) / _detD;
}



void Triangle::local_rhs_vector(double *loc_vec,
                                double(*rhs_func)(const Point &point, double t, const Parameters &par),
                                const std::vector<Point> &mesh_vertices, double time, const Parameters &param) const
{
  expect(fabs(_detD) > 1e-15,
         "An attempt to calculate a local matrix for a singular triangle (detD = " + d2s(_detD, true) + ")");

  switch (_dofs.size())
  {
  case n_dofs_first:
  {
    const double mat[][n_dofs_first] = { { 2., 1., 1. },
                                         { 1., 2., 1. },
                                         { 1., 1., 2. }
                                       };
    for (int i = 0; i < n_dofs_first; ++i)
    {
      loc_vec[i] = 0.;
      for (int j = 0; j < n_dofs_first; ++j)
        loc_vec[i] += mat[i][j] * rhs_func(mesh_vertices[_vertices[j]], time, param);
      loc_vec[i] *= fabs(_detD) / 24.;
    }
    return;
  }
  default:
    require(false, "Unknown fe order for generating local matrix");
  }
}
