#ifndef TESTING_H
#define TESTING_H

#include "config.h"
#include "fem/point.h"
#include "fem/triangle.h"
#include "fem/fine_mesh.h"
#include "auxiliary_testing_functions.h"
#include "fem/dof_handler.h"
#include "fem/csr_pattern.h"
#include "parameters.h"
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"
#include "fem/math_functions.h"
#include "acoustic2d.h"
#include <gtest/gtest.h>
#include <algorithm>
#include "analytic_functions.h"
#include "fem/function.h"


// =================================
//
// =================================
TEST(FineMesh, read_without_physical_and_partitions)
{
  fem::FineMesh fmesh;
  fmesh.read(TEST_DIR + "/test_mesh_0.msh",
             fem::Point(0, 0),
             fem::Point(1, 1));

  fem::Point points[] = { fem::Point(0, 0, 0),
                          fem::Point(1, 0, 0),
                          fem::Point(1, 1, 0),
                          fem::Point(0, 1, 0),
                          fem::Point(0.499999999998694, 0, 0),
                          fem::Point(1, 0.499999999998694, 0),
                          fem::Point(0.5000000000020591, 1, 0),
                          fem::Point(0, 0.5000000000020591, 0),
                          fem::Point(0.5, 0.5, 0),
                          fem::Point(0.2500000000010295, 0.7500000000010295, 0),
                          fem::Point(0.7500000000010295, 0.7499999999993469, 0),
                          fem::Point(0.7499999999993471, 0.249999999999347, 0),
                          fem::Point(0.2499999999994391, 0.2500000000006863, 0)
                        };

  fem::Triangle triangles[] ={ fem::Triangle(7, 4, 10),
                               fem::Triangle(1, 13, 8),
                               fem::Triangle(6, 3, 11),
                               fem::Triangle(2, 12, 5),
                               fem::Triangle(7, 10, 9),
                               fem::Triangle(8, 9, 10),
                               fem::Triangle(7, 9, 11),
                               fem::Triangle(8, 13, 9),
                               fem::Triangle(11, 9, 12),
                               fem::Triangle(5, 9, 13),
                               fem::Triangle(6, 11, 12),
                               fem::Triangle(5, 12, 9),
                               fem::Triangle(1, 5, 13),
                               fem::Triangle(6, 12, 2),
                               fem::Triangle(4, 8, 10),
                               fem::Triangle(7, 11, 3)
                             };

  compare_points(fmesh, 13, points);
  compare_triangles(fmesh, 16, triangles, points);
}



// =================================
//
// =================================
TEST(FineMesh, create_rectangular_grid)
{
}



// =================================
//
// =================================
TEST(FineMesh, read_with_physical_and_partitions)
{
  fem::FineMesh fmesh;
  fmesh.read(TEST_DIR + "/test_mesh_1.msh",
             fem::Point(0, 0),
             fem::Point(1, 1));

  // check the number of the mesh vertices
  EXPECT_EQ((int)fmesh.n_vertices(), 180);

  fem::Point points[] = { fem::Point(0, 0, 0), // first vertex
                     fem::Point(0.3207106781186477, 0.8500000000000002, 0) // last vertex
                   };

  // check the first vertex
  for (unsigned int j = 0; j < fem::Point::n_coord; ++j)
    EXPECT_DOUBLE_EQ(fmesh.vertex(0).coord(j), points[0].coord(j));

  // check the last vertex
  for (unsigned int j = 0; j < fem::Point::n_coord; ++j)
    EXPECT_DOUBLE_EQ(fmesh.vertex(fmesh.n_vertices() - 1).coord(j), points[1].coord(j));

  // check the number of the mesh triangles
  EXPECT_EQ((int)fmesh.n_triangles(), 318);

  // some vectors of "ghost cells" for triangles lying on partitions boundaries
  std::vector<unsigned int> gc1(1);
  gc1[0] = 2;

  std::vector<unsigned int> gc2(2);
  gc2[0] = 1;
  gc2[1] = 4;

  std::vector<unsigned int> gc3(1);
  gc3[0] = 5;

  std::vector<fem::Point> empty_points;

  fem::Triangle triangles[] ={ fem::Triangle(17, 89, 114, empty_points, 1, 1), // triangle number 1 (starting from 1) - first one
                               fem::Triangle(52, 115, 76, empty_points, 1, 5, gc1), // triangle number 9 (starting from 1)
                               fem::Triangle(99, 168, 121, empty_points, 1, 2, gc2), // triangle number 185 (starting from 1)
                               fem::Triangle(8, 47, 178, empty_points, 11, 4, gc3)  // triangle number 318 (starting from 1) - last one
                             };

  // triangle number 1
  for (unsigned int j = 0; j < fem::Triangle::n_vertices; ++j)
  {
    // check vertices numbers.
    // we do "+1" since initially the mesh vertices are numerated from 1. but in the program they are shifted to start from 0
    EXPECT_EQ(fmesh.triangle(0).vertex(j) + 1, triangles[0].vertex(j));
  }
  // check material id
  EXPECT_EQ(fmesh.triangle(0).material_id(), triangles[0].material_id());
  EXPECT_EQ((int)fmesh.triangle(0).material_id(), 1);
  // check partition id
  EXPECT_EQ(fmesh.triangle(0).partition_id() + 1, triangles[0].partition_id());
  EXPECT_EQ((int)fmesh.triangle(0).partition_id() + 1, 1);
  // check the number of ghost cells
  EXPECT_EQ(fmesh.triangle(0).n_ghost_cells(), triangles[0].n_ghost_cells());
  EXPECT_EQ((int)fmesh.triangle(0).n_ghost_cells(), 0);

  // triangle number 9
  for (unsigned int j = 0; j < fem::Triangle::n_vertices; ++j)
  {
    // check vertices numbers.
    // we do "+1" since initially the mesh vertices are numerated from 1. but in the program they are shifted to start from 0
    EXPECT_EQ(fmesh.triangle(8).vertex(j) + 1, triangles[1].vertex(j));
  }
  // check material id
  EXPECT_EQ(fmesh.triangle(8).material_id(), triangles[1].material_id());
  EXPECT_EQ((int)fmesh.triangle(8).material_id(), 1);
  // check partition id
  EXPECT_EQ(fmesh.triangle(8).partition_id() + 1, triangles[1].partition_id());
  EXPECT_EQ((int)fmesh.triangle(8).partition_id() + 1, 5);
  // check the number of ghost cells
  EXPECT_EQ(fmesh.triangle(8).n_ghost_cells(), triangles[1].n_ghost_cells());
  EXPECT_EQ(fmesh.triangle(8).n_ghost_cells(), gc1.size());
  // check ghost cells
  for (unsigned int g = 0; g < fmesh.triangle(8).n_ghost_cells(); ++g)
  {
    EXPECT_EQ(fmesh.triangle(8).ghost_cell(g) + 1, triangles[1].ghost_cell(g));
    EXPECT_EQ(fmesh.triangle(8).ghost_cell(g) + 1, gc1[g]);
  }

  // triangle number 185
  for (unsigned int j = 0; j < fem::Triangle::n_vertices; ++j)
  {
    // check vertices numbers.
    // we do "+1" since initially the mesh vertices are numerated from 1. but in the program they are shifted to start from 0
    EXPECT_EQ(fmesh.triangle(184).vertex(j) + 1, triangles[2].vertex(j));
  }
  // check material id
  EXPECT_EQ(fmesh.triangle(184).material_id(), triangles[2].material_id());
  EXPECT_EQ((int)fmesh.triangle(184).material_id(), 1);
  // check partition id
  EXPECT_EQ(fmesh.triangle(184).partition_id() + 1, triangles[2].partition_id());
  EXPECT_EQ((int)fmesh.triangle(184).partition_id() + 1, 2);
  // check the number of ghost cells
  EXPECT_EQ(fmesh.triangle(184).n_ghost_cells(), triangles[2].n_ghost_cells());
  EXPECT_EQ(fmesh.triangle(184).n_ghost_cells(), gc2.size());
  // check ghost cells
  for (unsigned int g = 0; g < fmesh.triangle(184).n_ghost_cells(); ++g)
  {
    EXPECT_EQ(fmesh.triangle(184).ghost_cell(g) + 1, triangles[2].ghost_cell(g));
    EXPECT_EQ(fmesh.triangle(184).ghost_cell(g) + 1, gc2[g]);
  }

  // triangle number 318
  for (unsigned int j = 0; j < fem::Triangle::n_vertices; ++j)
  {
    // check vertices numbers.
    // we do "+1" since initially the mesh vertices are numerated from 1. but in the program they are shifted to start from 0
    EXPECT_EQ(fmesh.triangle(317).vertex(j) + 1, triangles[3].vertex(j));
  }
  // check material id
  EXPECT_EQ(fmesh.triangle(317).material_id(), triangles[3].material_id());
  EXPECT_EQ((int)fmesh.triangle(317).material_id(), 11);
  // check partition id
  EXPECT_EQ(fmesh.triangle(317).partition_id() + 1, triangles[3].partition_id());
  EXPECT_EQ((int)fmesh.triangle(317).partition_id() + 1, 4);
  // check the number of ghost cells
  EXPECT_EQ(fmesh.triangle(317).n_ghost_cells(), triangles[3].n_ghost_cells());
  EXPECT_EQ(fmesh.triangle(317).n_ghost_cells(), gc3.size());
  // check ghost cells
  for (unsigned int g = 0; g < fmesh.triangle(317).n_ghost_cells(); ++g)
  {
    EXPECT_EQ(fmesh.triangle(317).ghost_cell(g) + 1, triangles[3].ghost_cell(g));
    EXPECT_EQ(fmesh.triangle(317).ghost_cell(g) + 1, gc3[g]);
  }

  gc1.clear();
  gc2.clear();
  gc3.clear();
}

#endif // TESTING_H




// =================================
//
// =================================
TEST(CSRPattern, check_pattern_and_other_stuff)
{
  fem::FineMesh fmesh;
  fmesh.read(TEST_DIR + "/test_small_10_4.msh",
             fem::Point(0, 0),
             fem::Point(4, 4));

  const unsigned int n_points = 10;
  fem::Point points[] ={ fem::Point(0, 0, 0),
                         fem::Point(2, 0, 0),
                         fem::Point(4, 0, 0),
                         fem::Point(1, 1, 0),
                         fem::Point(0, 2, 0),
                         fem::Point(2, 2, 0),
                         fem::Point(4, 2, 0),
                         fem::Point(0, 4, 0),
                         fem::Point(2, 4, 0),
                         fem::Point(4, 4, 0)
                       };
  compare_points(fmesh, n_points, points);

  fem::FiniteElement fe(1);

  fem::DoFHandler dof_handler(&fmesh);
  dof_handler.distribute_dofs(fe, fem::CG); // default parameters
  check_dof_handler(dof_handler, fmesh, 1); // 1 means first order basis functions

  fem::CSRPattern csr_pattern;
  csr_pattern.make_sparse_format(dof_handler, fem::CG);

  unsigned int row[] = { 0, 4, 10, 13, 18, 23, 30, 36, 40, 45, 48 };
  unsigned int col[] = { 0, 1, 3, 4,
                         0, 1, 2, 3, 5, 6,
                         1, 2, 6,
                         0, 1, 3, 4, 5,
                         0, 3, 4, 5, 7,
                         1, 3, 4, 5, 6, 7, 8,
                         1, 2, 5, 6, 8, 9,
                         4, 5, 7, 8,
                         5, 6, 7, 8, 9,
                         6, 8, 9
                       };

  check_csr_pattern(csr_pattern, n_points, &row[0], &col[0]);

}



// =================================
//
// =================================
TEST(CSRPattern, check_that_pattern_doesnt_fall)
{
  fem::FineMesh fmesh;
  fmesh.read(TEST_DIR + "/test_mesh_1.msh",
             fem::Point(0, 0),
             fem::Point(1, 1));

  fem::FiniteElement fe(1);

  fem::DoFHandler dof_handler(&fmesh);
  dof_handler.distribute_dofs(fe, fem::CG); // default parameters

  fem::CSRPattern csr_pattern;
  EXPECT_NO_THROW(csr_pattern.make_sparse_format(dof_handler, fem::CG)); // make a CSR pattern
}



// =================================
//
// =================================
TEST(Acoustic2D, check_fall_without_comp_environment)
{
  Parameters d_p;
  EXPECT_ANY_THROW(Acoustic2D task(&d_p)); // attempt to launch the solver with default parameters

  Parameters *de_p = new Parameters();
  EXPECT_ANY_THROW(Acoustic2D task(de_p)); // attempt to launch the solver with default parameters

//  Parameters *def_p;
//  EXPECT_ANY_THROW(Acoustic2D task(def_p)); // attempt to launch the solver with default parameters
}



// =================================
//
// =================================
const std::string test_mesh_0_files[] = { "test_mesh_0_05.msh",
                                          "test_mesh_0_025.msh",
                                          "test_mesh_0_0125.msh",
                                          "test_mesh_0_00625.msh",
                                          "test_mesh_0_003125.msh",
                                          "test_mesh_0_0015625.msh" };
const unsigned int N_times_tri_dense = 5; // the last mesh is too big to work with dense matrix
const unsigned int N_times_tri_sparse = 6;

TEST(EllipticAnalyticSolutionDense, AnalyticFunction_x_plus_y)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_tri_dense; ++i)
  {
    check_elliptic_solution_triangles(0, test_mesh_0_files[i], an_solution_1(), an_rhs_function_1(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}

TEST(EllipticAnalyticSolutionDense, AnalyticFunction_x_mult_y)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_tri_dense; ++i)
  {
    check_elliptic_solution_triangles(0, test_mesh_0_files[i], an_solution_2(), an_rhs_function_2(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}

TEST(EllipticAnalyticSolutionDense, AnalyticFunction_xx_plus_yy)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_tri_dense; ++i)
  {
    check_elliptic_solution_triangles(0, test_mesh_0_files[i], an_solution_3(), an_rhs_function_3(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}

TEST(EllipticAnalyticSolutionDense, AnalyticFunction_sinx_plus_siny)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_tri_dense; ++i)
  {
    check_elliptic_solution_triangles(0, test_mesh_0_files[i], an_solution_4(), an_rhs_function_4(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}

TEST(EllipticAnalyticSolutionDense, AnalyticFunction_expx)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_tri_dense; ++i)
  {
    check_elliptic_solution_triangles(0, test_mesh_0_files[i], an_solution_5(), an_rhs_function_5(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}



// =================================
//
// =================================
TEST(EllipticAnalyticSolutionSparse, AnalyticFunction_x_plus_y)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_tri_sparse; ++i)
  {
    check_elliptic_solution_triangles(1, test_mesh_0_files[i], an_solution_1(), an_rhs_function_1(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}

TEST(EllipticAnalyticSolutionSparse, AnalyticFunction_x_mult_y)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_tri_sparse; ++i)
  {
    check_elliptic_solution_triangles(1, test_mesh_0_files[i], an_solution_2(), an_rhs_function_2(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}

TEST(EllipticAnalyticSolutionSparse, AnalyticFunction_xx_plus_yy)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_tri_sparse; ++i)
  {
    check_elliptic_solution_triangles(1, test_mesh_0_files[i], an_solution_3(), an_rhs_function_3(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}

TEST(EllipticAnalyticSolutionSparse, AnalyticFunction_sinx_plus_siny)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_tri_sparse; ++i)
  {
    check_elliptic_solution_triangles(1, test_mesh_0_files[i], an_solution_4(), an_rhs_function_4(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}

TEST(EllipticAnalyticSolutionSparse, AnalyticFunction_expx)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_tri_sparse; ++i)
  {
    check_elliptic_solution_triangles(1, test_mesh_0_files[i], an_solution_5(), an_rhs_function_5(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}




// =================================
//
// =================================
const unsigned int N_times_rec_dense = 4; // the last mesh is too big to work with dense matrix
const unsigned int N_times_rec_sparse = 5;
const unsigned int n_beg = 4; // start number of cells in x- and y-directions

TEST(EllipticAnalyticSolutionRectanglesDense, AnalyticFunction_x_plus_y)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_rec_dense; ++i)
  {
    check_elliptic_solution_rectangles(0, pow(2, i)*n_beg, pow(2, i)*n_beg, an_solution_1(), an_rhs_function_1(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}

TEST(EllipticAnalyticSolutionRectanglesDense, AnalyticFunction_x_mult_y)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_rec_dense; ++i)
  {
    check_elliptic_solution_rectangles(0, pow(2, i)*n_beg, pow(2, i)*n_beg, an_solution_2(), an_rhs_function_2(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}

TEST(EllipticAnalyticSolutionRectanglesDense, AnalyticFunction_xx_plus_yy)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_rec_dense; ++i)
  {
    check_elliptic_solution_rectangles(0, pow(2, i)*n_beg, pow(2, i)*n_beg, an_solution_3(), an_rhs_function_3(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}

TEST(EllipticAnalyticSolutionRectanglesDense, AnalyticFunction_sinx_plus_siny)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_rec_dense; ++i)
  {
    check_elliptic_solution_rectangles(0, pow(2, i)*n_beg, pow(2, i)*n_beg, an_solution_4(), an_rhs_function_4(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}

TEST(EllipticAnalyticSolutionRectanglesDense, AnalyticFunction_expx)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_rec_dense; ++i)
  {
    check_elliptic_solution_rectangles(0, pow(2, i)*n_beg, pow(2, i)*n_beg, an_solution_5(), an_rhs_function_5(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}



// =================================
//
// =================================
TEST(EllipticAnalyticSolutionRectanglesSparse, AnalyticFunction_x_plus_y)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_rec_sparse; ++i)
  {
    check_elliptic_solution_rectangles(1, pow(2, i)*n_beg, pow(2, i)*n_beg, an_solution_1(), an_rhs_function_1(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}

TEST(EllipticAnalyticSolutionRectanglesSparse, AnalyticFunction_x_mult_y)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_rec_sparse; ++i)
  {
    check_elliptic_solution_rectangles(1, pow(2, i)*n_beg, pow(2, i)*n_beg, an_solution_2(), an_rhs_function_2(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}

TEST(EllipticAnalyticSolutionRectanglesSparse, AnalyticFunction_xx_plus_yy)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_rec_sparse; ++i)
  {
    check_elliptic_solution_rectangles(1, pow(2, i)*n_beg, pow(2, i)*n_beg, an_solution_3(), an_rhs_function_3(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}

TEST(EllipticAnalyticSolutionRectanglesSparse, AnalyticFunction_sinx_plus_siny)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_rec_sparse; ++i)
  {
    check_elliptic_solution_rectangles(1, pow(2, i)*n_beg, pow(2, i)*n_beg, an_solution_4(), an_rhs_function_4(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}

TEST(EllipticAnalyticSolutionRectanglesSparse, AnalyticFunction_expx)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;
  for (unsigned int i = 0; i < N_times_rec_sparse; ++i)
  {
    check_elliptic_solution_rectangles(1, pow(2, i)*n_beg, pow(2, i)*n_beg, an_solution_5(), an_rhs_function_5(), cur_error, cur_n_cells, prev_error, prev_n_cells);
    prev_error = cur_error;
    prev_n_cells = cur_n_cells;
  }
}



// =================================
//
// =================================
TEST(EllipticAnalyticSolutionRectanglesDiffDomainsAndMeshes, AnalyticFunction_expx)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;

  // sparse solvers
  check_elliptic_solution_rectangles(1, 50, 50, an_solution_4(), an_rhs_function_4(),
                                     cur_error, cur_n_cells, prev_error, prev_n_cells,
                                     0, 10, 0, 10);
  check_elliptic_solution_rectangles(1, 51, 51, an_solution_4(), an_rhs_function_4(),
                                     cur_error, cur_n_cells, prev_error, prev_n_cells,
                                     0, 10, 0, 10);
  check_elliptic_solution_rectangles(1, 53, 53, an_solution_4(), an_rhs_function_4(),
                                     cur_error, cur_n_cells, prev_error, prev_n_cells,
                                     0, 10, 0, 10);
//  check_elliptic_solution_rectangles(1, 10, 11, an_solution_4(), an_rhs_function_4(),
//                                     cur_error, cur_n_cells, prev_error, prev_n_cells,
//                                     0, 1000, 0, 1000);
//  check_elliptic_solution_rectangles(1, 11, 10, an_solution_4(), an_rhs_function_4(),
//                                     cur_error, cur_n_cells, prev_error, prev_n_cells,
//                                     0, 1000, 0, 1000);
//  check_elliptic_solution_rectangles(1, 13, 13, an_solution_4(), an_rhs_function_4(),
//                                     cur_error, cur_n_cells, prev_error, prev_n_cells,
//                                     0, 1000, 0, 1000);
//  check_elliptic_solution_rectangles(1, 131, 131, an_solution_4(), an_rhs_function_4(),
//                                     cur_error, cur_n_cells, prev_error, prev_n_cells,
//                                     1000, 1500, 1000, 1500);
//  check_elliptic_solution_rectangles(1, 151, 131, an_solution_4(), an_rhs_function_4(),
//                                     cur_error, cur_n_cells, prev_error, prev_n_cells,
//                                     1000, 1500, 1000, 1500);
//  check_elliptic_solution_rectangles(1, 131, 151, an_solution_4(), an_rhs_function_4(),
//                                     cur_error, cur_n_cells, prev_error, prev_n_cells,
//                                     1000, 1500, 1000, 1500);
}



// =================================
//
// =================================
TEST(EllipticAnalyticSolutionRectanglesConstants, AnalyticFunction_constant_1)
{
  double cur_error, prev_error = -1;
  int cur_n_cells, prev_n_cells = 0;

  // sparse solvers
  check_elliptic_solution_rectangles(1, 2, 2, fem::ConstantFunction(1), fem::ConstantFunction(0),
                                     cur_error, cur_n_cells, prev_error, prev_n_cells,
                                     0, 10, 0, 10);
}
