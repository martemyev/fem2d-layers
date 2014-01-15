#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "mesh_element.h"
#include "point.h"

/**
 * Triangle - 2-dimensional simplex.
 * The simplest shape in 2D.
 * It's an element of mesh,
 * therefore it inherits the most part of
 * functionality from MeshElement class.
 */
class Triangle : public MeshElement
{
public:
            /**
             * The number of vertices of triangle
             */
  static const unsigned int n_vertices = 3;

            /**
             * The number of edges of triangle
             */
  static const unsigned int n_edges = 3;

            /**
             * Triangle is 2D shape,
             * so it's a face itself (for tetrahedron)
             */
  static const unsigned int n_faces = 1;

            /**
             * In Gmsh triangle is defined by number 2
             */
  static const unsigned int gmsh_el_type = 2;

            /**
             * In VTK format triangle is defined by number 5
             */
  static const unsigned int vtk_el_type = 5;

            /**
             * The number of degrees of freedom on the triangle
             * in case of first order basis functions coincides
             * with the number of vertices
             */
  static const unsigned int n_dofs_first = n_vertices;

            /**
             * Default constructor
             */
  Triangle();

            /**
             * Destructor
             */
  ~Triangle();

            /**
             * Constructor with parameters
             * @param ver - triangle vertices
             * @param mat_id - material ID
             */
  Triangle(const std::vector<unsigned int> &ver,
           const std::vector<Point> &mesh_vertices = std::vector<Point>(),
           const unsigned int mat_id = 0,
           const unsigned int part_id = 0,
           const std::vector<unsigned int> &ghost_cells = std::vector<unsigned int>());

            /**
             * Constructor with parameters
             * @param v1 - first vertex
             * @param v2 - second vertex
             * @param v3 - third vertex
             * @param mat_id - material ID
             */
  Triangle(const unsigned int v1,
           const unsigned int v2,
           const unsigned int v3,
           const std::vector<Point> &mesh_vertices = std::vector<Point>(),
           const unsigned int mat_id = 0,
           const unsigned int part_id = 0,
           const std::vector<unsigned int> &ghost_cells = std::vector<unsigned int>());

            /**
             * Copy constructor
             */
  Triangle(const Triangle &tri);

            /**
             * Copy assignment operator
             */
  Triangle& operator =(const Triangle &tri);

            /**
             * Get the number of ghost cells
             */
  unsigned int n_ghost_cells() const;

            /**
             * Get the number of 'num'-th ghost cell
             * @param num - the serial number of ghost cell fomr the vector
             */
  unsigned int ghost_cell(unsigned int num) const;

            /**
             * Get the number of degrees of freedom associated with this triangle
             */
  unsigned int n_dofs() const;

            /**
             * Get the number of degree of freedom
             * @param number - the serial number of the dof for this triangle
             */
  unsigned int dof(unsigned int number) const;

            /**
             * Set the number of degrees of freedom on the triangle.
             * Actually - just to allocate the memory for _dofs vector
             * @param number - the number of degrees of freedom associated with the triangle
             */
  void n_dofs(unsigned int number);

            /**
             * Set the value of a degree of freedom
             * @param number - the serial number of the dof
             * @param value - the value of the dof
             */
  void dof(unsigned int number, unsigned int value);

            /**
             * Generate the local mass matrix for the triangle
             * @param loc_mat - output data, generated local matrix
             */
  void local_mass_matrix(double **loc_mat) const;

            /**
             * Generate the local stiffness matrix for the triangle
             * @param loc_mat - output data, generated local matrix
             */
  void local_stiffness_matrix(double **loc_mat) const;

            /**
             * Generate the local rhs vector
             * @param loc_vec - output vector
             * @param rhs_func - pointer to an rhs function
             */
  void local_rhs_vector(double *loc_vec, double(*rhs_func)(const Point &point, double t), const std::vector<Point> &mesh_vertices, double time = 0.) const;


private:
            /**
             * The number of partitions which the triangles is connected with.
             */
  std::vector<unsigned int> _ghost_cells;

            /**
             * Degrees of freedom defined on this triangle
             */
  std::vector<unsigned int> _dofs;

            /**
             * Determinant of the special matrix of coordinates of the triangles vertices
             */
  double _detD;

            /**
             * Coefficients _A of the barycentric functions on the triangle:
             * L[i] = _A[i]*x + _B[i]*y + _C[i]
             */
  double _A[n_vertices];

            /**
             * Coefficients _B of the barycentric functions on the triangle:
             * L[i] = _A[i]*x + _B[i]*y + _C[i]
             */
  double _B[n_vertices];

            /**
             * Coefficients _C of the barycentric functions on the triangle:
             * L[i] = _A[i]*x + _B[i]*y + _C[i]
             */
  double _C[n_vertices];

            /**
             * Calculate all coefficients of barycentric functions and detD
             * @param mesh_vertices - the vertices of the mesh
             */
  void calculate_barycentric(const std::vector<Point> &mesh_vertices);
};




#endif // TRIANGLE_H
