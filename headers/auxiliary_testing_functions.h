#ifndef AUXILARY_TESTING_FUNCTIONS_H
#define AUXILARY_TESTING_FUNCTIONS_H

#include "fem/fine_mesh.h"
#include "fem/point.h"
#include "fem/triangle.h"
#include "fem/dof_handler.h"
#include "fem/csr_pattern.h"
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"
#include "fem/math_functions.h"
#include "fem/auxiliary_functions.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <iostream>


// solver parameters
const double ksp_rtol = 1e-12;
const double ksp_atol = 1e-30;
const double ksp_dtol = 1e+5;
const int ksp_maxits  = 10000;


void compare_points(const FineMesh &fmesh, unsigned int n_points, Point points[])
{
  // check the number of vertices
  EXPECT_EQ(fmesh.n_vertices(), n_points);

  // check the coordinates
  for (int i = 0; i < fmesh.n_vertices(); ++i)
    for (int j = 0; j < Point::n_coord; ++j)
      EXPECT_DOUBLE_EQ(fmesh.vertex(i).coord(j), points[i].coord(j));
}



void compare_triangles(const FineMesh &fmesh, unsigned int n_triangles, Triangle triangles[], Point points[])
{
  // check the number of triangles
  EXPECT_EQ(fmesh.n_triangles(), n_triangles);

  // check triangles themselves
  for (int i = 0; i < fmesh.n_triangles(); ++i)
  {
    for (int j = 0; j < Triangle::n_vertices; ++j)
    {
      // check vertices numbers.
      // we do "+1" since initially the mesh vertices are numerated from 1. but in the program they are shifted to start from 0
      EXPECT_EQ(fmesh.triangle(i).vertex(j) + 1, triangles[i].vertex(j));

      // check vertices coordinates (once again)
      for (int k = 0; k < Point::n_coord; ++k)
        EXPECT_DOUBLE_EQ(fmesh.vertex(fmesh.triangle(i).vertex(j)).coord(k), points[triangles[i].vertex(j) - 1].coord(k));
    }
  }
}



void check_dof_handler(const DoFHandler &dof_handler, const FineMesh &fmesh, unsigned int fe_order)
{
  if (fe_order == 1)
  {
    EXPECT_EQ(dof_handler.n_dofs(), fmesh.n_vertices());
    for (int i = 0; i < dof_handler.n_dofs(); ++i)
      for (int j = 0; j < Point::n_coord; ++j)
        EXPECT_DOUBLE_EQ(dof_handler.dof(i).coord(j), fmesh.vertex(i).coord(j));
  }
}



void check_csr_pattern(const CSRPattern &csr_pattern, unsigned int n_points, unsigned int *row, unsigned int *col)
{
  EXPECT_EQ(csr_pattern.order(), n_points);

  for (int i = 0; i < n_points + 1; ++i)
    EXPECT_EQ(csr_pattern.row(i), row[i]);

  for (int i = 0; i < csr_pattern.row(n_points); ++i)
    EXPECT_EQ(csr_pattern.col(i), col[i]);
}



double check_elliptic_solution_dense(const std::string &meshfile,
                                     double(*an_solution)(const Point &p, double t),
                                     double(*rhs_function)(const Point &p, double t, const Parameters &par),
                                     int prev_n_triangles = 0,
                                     double prev_rel_error = -1)
{
  Parameters default_param;

  FineMesh fmesh(&default_param);
  fmesh.read(TEST_DIR + "/" + meshfile);

  DoFHandler dof_handler(&fmesh);
  dof_handler.distribute_dofs(default_param, CG);

  // create vectors
  Vec system_rhs; // right hand side vector
  Vec solution; // numerical solution
  Vec exact_solution; // analytic solution

  // allocate memory
  VecCreateSeq(PETSC_COMM_SELF, fmesh.n_vertices(), &system_rhs);
  VecDuplicate(system_rhs, &solution);
  VecDuplicate(system_rhs, &exact_solution);

  // create PETSc matrix
  Mat dense_stiff_mat;
  MatCreateDense(PETSC_COMM_WORLD, dof_handler.n_dofs(), dof_handler.n_dofs(), dof_handler.n_dofs(), dof_handler.n_dofs(), NULL, &dense_stiff_mat);

  // allocate the memory for local matrices and vectors
  double **local_stiff_mat = new double*[Triangle::n_dofs_first];
  double *local_rhs_vec = new double[Triangle::n_dofs_first];
  for (int i = 0; i < Triangle::n_dofs_first; ++i)
    local_stiff_mat[i] = new double[Triangle::n_dofs_first];

  const double time = 0.; // because it's elliptic problem

  // assemble the matrix
  for (int cell = 0; cell < fmesh.n_triangles(); ++cell)
  {
    Triangle triangle = fmesh.triangle(cell);
    triangle.local_stiffness_matrix(1, local_stiff_mat); // 1 is the coefficient a
    triangle.local_rhs_vector(rhs_function, fmesh.vertices(), time, default_param, local_rhs_vec);

    for (int i = 0; i < triangle.n_dofs(); ++i)
    {
      const int dof_i = triangle.dof(i);
      VecSetValue(system_rhs, dof_i, local_rhs_vec[i], ADD_VALUES);
      for (int j = 0; j < triangle.n_dofs(); ++j)
      {
        const int dof_j = triangle.dof(j);
        MatSetValue(dense_stiff_mat, dof_i, dof_j, local_stiff_mat[i][j], ADD_VALUES);
      }
    }
  }

  // free the memory
  for (int i = 0; i < Triangle::n_dofs_first; ++i)
    delete[] local_stiff_mat[i];
  delete[] local_stiff_mat;
  delete[] local_rhs_vec;

  MatAssemblyBegin(dense_stiff_mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(dense_stiff_mat, MAT_FINAL_ASSEMBLY);

  // Dirichlet bc

  // find max and min values to know the limits of computational domain.
  // this approach works only if you have a rectangular domain
  double xmin = fmesh.vertex(0).coord(0);
  double xmax = fmesh.vertex(0).coord(0);
  double ymin = fmesh.vertex(0).coord(1);
  double ymax = fmesh.vertex(0).coord(1);
  for (int i = 1; i < fmesh.n_vertices(); ++i)
  {
    const double x = fmesh.vertex(i).coord(0);
    const double y = fmesh.vertex(i).coord(1);
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
  }

  // find nodes that lie on these limits
  std::vector<int> b_nodes;
  for (int i = 0; i < fmesh.n_vertices(); ++i)
  {
    const double x = fmesh.vertex(i).coord(0);
    const double y = fmesh.vertex(i).coord(1);
    if ((fabs(x - xmin) < 1e-14 || fabs(x - xmax) < 1e-14 ||
         fabs(y - ymin) < 1e-14 || fabs(y - ymax) < 1e-14) &&
        find(b_nodes.begin(), b_nodes.end(), i) == b_nodes.end())
      b_nodes.push_back(i);
  }

  // impose Dirichlet boundary condition
  // with ones on diagonal
  MatZeroRows(dense_stiff_mat, b_nodes.size(), &b_nodes[0], 1., solution, system_rhs); // change the matrix
  for (int i = 0; i < b_nodes.size(); ++i)
    VecSetValue(system_rhs, b_nodes[i], an_solution(fmesh.vertex(b_nodes[i]), time), INSERT_VALUES); // change the rhs vector

  // solve the SLAE
  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, dense_stiff_mat, dense_stiff_mat, SAME_PRECONDITIONER);
  KSPSetTolerances(ksp, ksp_rtol, ksp_atol, ksp_dtol, ksp_maxits);
  KSPSolve(ksp, system_rhs, solution);

  // check solution
  for (int i = 0; i < fmesh.n_vertices(); ++i)
  {
    Point vert = fmesh.vertex(i);
    VecSetValue(exact_solution, i, an_solution(vert, time), INSERT_VALUES);
  }

  const int n_triangles = fmesh.n_triangles(); // the number of triangles of fine mesh
  const double rel_err = rel_error(solution, exact_solution); // relative error
  std::cout << "relative error = " << d2s(rel_err, true)
            << " error reduction = " << (fabs(prev_rel_error + 1) < 1e-12 ? "-" : d2s(prev_rel_error / rel_err))
            << " mesh size reduction = " << (prev_n_triangles == 0 ? "-" : d2s(n_triangles / prev_n_triangles))
            << std::endl;

  KSPDestroy(&ksp);

  MatDestroy(&dense_stiff_mat);

  VecDestroy(&system_rhs);
  VecDestroy(&solution);
  VecDestroy(&exact_solution);

  return rel_err;
}



double check_elliptic_solution_sparse(const std::string &meshfile,
                                      double(*an_solution)(const Point &p, double t),
                                      double(*rhs_function)(const Point &p, double t, const Parameters &par),
                                      int prev_n_triangles = 0,
                                      double prev_rel_error = -1)
{
  Parameters default_param;
  FineMesh fmesh(&default_param);
  fmesh.read(TEST_DIR + "/" + meshfile);
  DoFHandler dof_handler(&fmesh);
  dof_handler.distribute_dofs(default_param, CG);
  CSRPattern csr_pattern;
  csr_pattern.make_sparse_format(dof_handler, CG);

  // create vectors
  Vec system_rhs; // right hand side vector
  Vec solution; // numerical solution
  Vec exact_solution; // analytic solution

  // allocate memory
  VecCreateSeq(PETSC_COMM_SELF, fmesh.n_vertices(), &system_rhs);
  VecDuplicate(system_rhs, &solution);
  VecDuplicate(system_rhs, &exact_solution);

  // create PETSc matrix
  Mat sparse_stiff_mat;
  MatCreateSeqAIJ(PETSC_COMM_WORLD, csr_pattern.order(), csr_pattern.order(), 0, csr_pattern.nnz(), &sparse_stiff_mat);

  // allocate the memory for local matrices and vectors
  double **local_stiff_mat = new double*[Triangle::n_dofs_first];
  double *local_rhs_vec = new double[Triangle::n_dofs_first];
  for (int i = 0; i < Triangle::n_dofs_first; ++i)
    local_stiff_mat[i] = new double[Triangle::n_dofs_first];

  const double time = 0.; // because it's elliptic problem

  // assemble the matrix
  for (int cell = 0; cell < fmesh.n_triangles(); ++cell)
  {
    Triangle triangle = fmesh.triangle(cell);
    triangle.local_stiffness_matrix(1, local_stiff_mat); // 1 is the coefficient a
    triangle.local_rhs_vector(rhs_function, fmesh.vertices(), time, default_param, local_rhs_vec);

    for (int i = 0; i < triangle.n_dofs(); ++i)
    {
      const int dof_i = triangle.dof(i);
      VecSetValue(system_rhs, dof_i, local_rhs_vec[i], ADD_VALUES);
      for (int j = 0; j < triangle.n_dofs(); ++j)
      {
        const int dof_j = triangle.dof(j);
        MatSetValue(sparse_stiff_mat, dof_i, dof_j, local_stiff_mat[i][j], ADD_VALUES);
      }
    }
  }

  // free the memory
  for (int i = 0; i < Triangle::n_dofs_first; ++i)
    delete[] local_stiff_mat[i];
  delete[] local_stiff_mat;
  delete[] local_rhs_vec;

  MatAssemblyBegin(sparse_stiff_mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(sparse_stiff_mat, MAT_FINAL_ASSEMBLY);

  // Dirichlet bc

  // find max and min values to know the limits of computational domain.
  // this approach works only if you have a rectangular domain
  double xmin = fmesh.vertex(0).coord(0);
  double xmax = fmesh.vertex(0).coord(0);
  double ymin = fmesh.vertex(0).coord(1);
  double ymax = fmesh.vertex(0).coord(1);
  for (int i = 1; i < fmesh.n_vertices(); ++i)
  {
    const double x = fmesh.vertex(i).coord(0);
    const double y = fmesh.vertex(i).coord(1);
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
  }

  // find nodes that lie on these limits
  std::vector<int> b_nodes;
  for (int i = 0; i < fmesh.n_vertices(); ++i)
  {
    const double x = fmesh.vertex(i).coord(0);
    const double y = fmesh.vertex(i).coord(1);
    if ((fabs(x - xmin) < 1e-14 || fabs(x - xmax) < 1e-14 ||
         fabs(y - ymin) < 1e-14 || fabs(y - ymax) < 1e-14) &&
        find(b_nodes.begin(), b_nodes.end(), i) == b_nodes.end())
      b_nodes.push_back(i);
  }

  // impose Dirichlet boundary condition
  // with ones on diagonal
  MatZeroRows(sparse_stiff_mat, b_nodes.size(), &b_nodes[0], 1., solution, system_rhs); // change the matrix
  for (int i = 0; i < b_nodes.size(); ++i)
    VecSetValue(system_rhs, b_nodes[i], an_solution(fmesh.vertex(b_nodes[i]), time), INSERT_VALUES); // change the rhs vector

  // solve the SLAE
  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, sparse_stiff_mat, sparse_stiff_mat, SAME_PRECONDITIONER);
  KSPSetTolerances(ksp, ksp_rtol, ksp_atol, ksp_dtol, ksp_maxits);
  KSPSolve(ksp, system_rhs, solution);

  // check solution
  for (int i = 0; i < fmesh.n_vertices(); ++i)
  {
    Point vert = fmesh.vertex(i);
    VecSetValue(exact_solution, i, an_solution(vert, time), INSERT_VALUES);
  }

  const int n_triangles = fmesh.n_triangles(); // the number of triangles of fine mesh
  const double rel_err = rel_error(solution, exact_solution); // relative error
  std::cout << "relative error = " << d2s(rel_err, true)
            << " error reduction = " << (fabs(prev_rel_error + 1) < 1e-12 ? "-" : d2s(prev_rel_error / rel_err))
            << " mesh size reduction = " << (prev_n_triangles == 0 ? "-" : d2s(n_triangles / prev_n_triangles))
            << std::endl;

  KSPDestroy(&ksp);

  MatDestroy(&sparse_stiff_mat);

  VecDestroy(&system_rhs);
  VecDestroy(&solution);
  VecDestroy(&exact_solution);

  return rel_err;
}



#endif // AUXILARY_TESTING_FUNCTIONS_H
