#include "acoustic2d.h"
#include "parameters.h"
#include "petscksp.h"
#include "auxiliary_functions.h"
#include "analytic_functions.h"
#include "result.h"
#include "math_functions.h"
#include <iostream>
#include <algorithm>


Acoustic2D::Acoustic2D(Parameters *param)
  : _param(param)
{ }


Acoustic2D::~Acoustic2D()
{
  //VecDestroy(&_global_rhs);
  //MatDestroy(&_global_mass_mat);
  //MatDestroy(&_global_stiff_mat);
}



void Acoustic2D::solve()
{
  require(_param->FE_ORDER == 1, "This fe order is not implemented (" + d2s(_param->FE_ORDER) + ")");

  // read the fine triangular mesh from the file
  _fmesh.read(_param->MESH_FILE);
#if defined(DEBUG)
  std::cout << "n_nodes = " << _fmesh.n_vertices() << std::endl;
#endif

  DoFHandler dof_handler(&_fmesh, *_param);

  // create sparse format based on the distribution of degrees of freedom.
  // since we use first order basis functions, and then
  // all dofs are associated with the mesh vertices,
  // sparse format is based on connectivity of the mesh vertices
  CSRPattern csr_pattern(dof_handler);

  expect(csr_pattern.order() == dof_handler.n_dofs(), "Error");
#if defined(DEBUG)
  std::cout << "n_dofs = " << dof_handler.n_dofs() << std::endl;
#endif

  // allocate memory
  VecCreateSeq(PETSC_COMM_SELF, csr_pattern.order(), &_global_rhs);
  MatCreateSeqAIJ(PETSC_COMM_WORLD, csr_pattern.order(), csr_pattern.order(), 0, csr_pattern.nnz(), &_global_mass_mat);
  MatCreateSeqAIJ(PETSC_COMM_WORLD, csr_pattern.order(), csr_pattern.order(), 0, csr_pattern.nnz(), &_global_stiff_mat);

  // allocate the memory for local matrices and vectors
  double **local_mass_mat = new double*[Triangle::n_dofs_first];
  double **local_stiff_mat = new double*[Triangle::n_dofs_first];
  for (int i = 0; i < Triangle::n_dofs_first; ++i)
  {
    local_mass_mat[i] = new double[Triangle::n_dofs_first];
    local_stiff_mat[i] = new double[Triangle::n_dofs_first];
  }

  // assemble the matrices and the rhs vector
  for (int cell = 0; cell < _fmesh.n_triangles(); ++cell)
  {
    Triangle triangle = _fmesh.triangle(cell);
    triangle.local_mass_matrix(local_mass_mat);
    triangle.local_stiffness_matrix(local_stiff_mat);

    for (int i = 0; i < triangle.n_dofs(); ++i)
    {
      const int dof_i = triangle.dof(i);
      //VecSetValue(_global_rhs, dof_i, local_rhs_vec[i], ADD_VALUES);
      for (int j = 0; j < triangle.n_dofs(); ++j)
      {
        const int dof_j = triangle.dof(j);
        MatSetValue(_global_mass_mat, dof_i, dof_j, local_mass_mat[i][j], ADD_VALUES);
        MatSetValue(_global_stiff_mat, dof_i, dof_j, local_stiff_mat[i][j], ADD_VALUES);
      }
    }
  }

  // free the memory
  for (int i = 0; i < Triangle::n_dofs_first; ++i)
    delete[] local_stiff_mat[i];
  delete[] local_stiff_mat;

  MatAssemblyBegin(_global_mass_mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(_global_mass_mat, MAT_FINAL_ASSEMBLY);

  MatAssemblyBegin(_global_stiff_mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(_global_stiff_mat, MAT_FINAL_ASSEMBLY);

  if (_param->TIME_SCHEME == EXPLICIT)
    solve_explicit(dof_handler, csr_pattern);
  else if(_param->TIME_SCHEME == CRANK_NICOLSON)
    solve_crank_nicolson(dof_handler, csr_pattern);
  else
    require(false, "Unknown time discretization scheme");


}



void Acoustic2D::solve_crank_nicolson(const DoFHandler &dof_handler, const CSRPattern &csr_pattern)
{ }



void Acoustic2D::solve_explicit(const DoFHandler &dof_handler, const CSRPattern &csr_pattern)
{
  // create vectors
  Vec solution; // numerical solution on the current (n-th) time step
  Vec solution_1; // numerical solution on the previous ((n-1)-th) time step
  Vec solution_2; // numerical solution on the preprevious ((n-2)-th) time step
  Vec exact_solution; // analytic solution

  VecDuplicate(_global_rhs, &solution);
  VecDuplicate(_global_rhs, &solution_1);
  VecDuplicate(_global_rhs, &solution_2);
  VecDuplicate(_global_rhs, &exact_solution);

  require(_param->FE_ORDER == 1, "This fe order is not implemented (" + d2s(_param->FE_ORDER) + ")");

  const double dt = _param->TIME_STEP;

  // fill vectors with solution on the 0-th and 1-st time steps
  for (int d = 0; d < dof_handler.n_dofs(); ++d)
  {
    VecSetValue(solution_2, d, an_init_solution(dof_handler.dof(d), _param->TIME_BEG), INSERT_VALUES);
    VecSetValue(solution_1, d, an_init_solution(dof_handler.dof(d), _param->TIME_BEG + dt), INSERT_VALUES);
  }

  // make a SLAE rhs vector
  Vec system_rhs, temp;
  VecCreateSeq(PETSC_COMM_SELF, csr_pattern.order(), &system_rhs);
  VecDuplicate(system_rhs, &temp);

  // boundary nodes
  std::vector<int> b_nodes;
  find_bound_nodes(b_nodes);

  double *local_rhs_vec = new double[Triangle::n_dofs_first];

  require(_param->N_TIME_STEPS > 2, "There is no time steps to perform: n_time_steps = " + d2s(_param->N_TIME_STEPS));
  for (int time_step = 2; time_step < _param->N_TIME_STEPS; ++time_step)
  {
    const double time = _param->TIME_BEG + time_step * dt; // current time
    VecSet(system_rhs, 0.); // zeroing the system rhs vector
    VecSet(solution, 0.); // zeroing the solution vector

    // assemble some parts of system rhs vector
    for (int cell = 0; cell < _fmesh.n_triangles(); ++cell)
    {
      Triangle triangle = _fmesh.triangle(cell);
      triangle.local_rhs_vector(local_rhs_vec, rhs_function, _fmesh.vertices(), time - dt); // rhs function on the previous time step
      for (int i = 0; i < triangle.n_dofs(); ++i)
      {
        const int dof_i = triangle.dof(i);
        VecSetValue(system_rhs, dof_i, local_rhs_vec[i], ADD_VALUES);
      }
    } // rhs part assembling

    MatMult(_global_stiff_mat, solution_1, temp);
    VecAXPBY(system_rhs, -dt*dt, dt*dt, temp);

    MatMult(_global_mass_mat, solution_2, temp);
    VecAXPY(system_rhs, -1., temp);

    MatMult(_global_mass_mat, solution_1, temp);
    VecAXPY(system_rhs, 2., temp);

    // system matrix equal to global mass matrix
    Mat system_mat;
    MatConvert(_global_mass_mat, MATSAME, MAT_INITIAL_MATRIX, &system_mat); // allocate memory and copy values

    // impose Dirichlet boundary condition
    // with ones on diagonal
    MatZeroRows(system_mat, b_nodes.size(), &b_nodes[0], 1., solution, system_rhs); // change the matrix
    for (int i = 0; i < b_nodes.size(); ++i)
      VecSetValue(system_rhs, b_nodes[i], an_bound_solution(_fmesh.vertex(b_nodes[i]), time), INSERT_VALUES); // change the rhs vector

    // solve the SLAE
    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, system_mat, system_mat, SAME_PRECONDITIONER);
    KSPSetTolerances(ksp, 1e-12, 1e-30, 1e+5, 10000);
    KSPSolve(ksp, system_rhs, solution);
    KSPDestroy(&ksp);

    MatDestroy(&system_mat);

    // check solution
//    for (int i = 0; i < _fmesh.n_vertices(); ++i)
//    {
//      Point vert = _fmesh.vertex(i);
//      VecSetValue(exact_solution, i, an_solution(vert, time), INSERT_VALUES);
//    }
//    std::cout << "time step = " << time_step << " time = " << time << " relative error = " << rel_error(solution, exact_solution) << std::endl;

    Result res(&dof_handler);
    std::string fname = "res/results-" + d2s(time_step) + ".vtu";
    res.write_vtu(fname, solution); //, exact_solution);

    // reassign the solutions on the previuos time steps
    VecCopy(solution_1, solution_2);
    VecCopy(solution,   solution_1);

  } // time loop

  delete[] local_rhs_vec;

  VecDestroy(&solution);
  VecDestroy(&solution_1);
  VecDestroy(&solution_2);
  VecDestroy(&exact_solution);
  VecDestroy(&system_rhs);
  VecDestroy(&temp);
}



//// create vectors
//Vec system_rhs; // right hand side vector
//Vec solution;

//VecCreateSeq(PETSC_COMM_SELF, csr_pattern.order(), &system_rhs);
////VecSet(system_rhs, 0.);

//VecCreateSeq(PETSC_COMM_SELF, csr_pattern.order(), &solution);
//VecSet(solution, 0.);

//// create PETSc matrix
//Mat dense_mass_mat;
//Mat dense_stiff_mat;
//MatCreateDense(PETSC_COMM_WORLD, csr_pattern.order(), csr_pattern.order(), csr_pattern.order(), csr_pattern.order(), NULL, &dense_mass_mat);
//MatCreateDense(PETSC_COMM_WORLD, csr_pattern.order(), csr_pattern.order(), csr_pattern.order(), csr_pattern.order(), NULL, &dense_stiff_mat);


//// assemble the matrix
//double **local_mass_mat = new double*[Triangle::n_dofs_first];
//double **local_stiff_mat = new double*[Triangle::n_dofs_first];
//double *local_rhs_vec = new double[Triangle::n_dofs_first];
//for (int i = 0; i < Triangle::n_dofs_first; ++i)
//{
//  local_mass_mat[i] = new double[Triangle::n_dofs_first];
//  local_stiff_mat[i] = new double[Triangle::n_dofs_first];
//}

//for (int cell = 0; cell < _fmesh.n_triangles(); ++cell)
//{
//  Triangle triangle = _fmesh.triangle(cell);
//  triangle.local_mass_matrix(local_mass_mat);
//  triangle.local_stiffness_matrix(local_stiff_mat);
//  triangle.local_rhs_vector(local_rhs_vec, rhs_function, _fmesh.vertices());

//  for (int i = 0; i < triangle.n_dofs(); ++i)
//  {
//    const int dof_i = triangle.dof(i);
//    VecSetValue(system_rhs, dof_i, local_rhs_vec[i], ADD_VALUES);
//    for (int j = 0; j < triangle.n_dofs(); ++j)
//    {
//      const int dof_j = triangle.dof(j);
//      MatSetValue(dense_mass_mat, dof_i, dof_j, local_mass_mat[i][j], ADD_VALUES);
//      MatSetValue(dense_stiff_mat, dof_i, dof_j, local_stiff_mat[i][j], ADD_VALUES);
//    }
//  }
//}

////VecView(system_rhs, PETSC_VIEWER_STDOUT_SELF);

//MatAssemblyBegin(dense_mass_mat, MAT_FINAL_ASSEMBLY);
//MatAssemblyEnd(dense_mass_mat, MAT_FINAL_ASSEMBLY);

//MatAssemblyBegin(dense_stiff_mat, MAT_FINAL_ASSEMBLY);
//MatAssemblyEnd(dense_stiff_mat, MAT_FINAL_ASSEMBLY);

//// Dirichlet bc

//// find max and min values to know the limits of computational domain.
//// this approach works only if you have a rectangular domain
//double xmin = _fmesh.vertex(0).coord(0);
//double xmax = _fmesh.vertex(0).coord(0);
//double ymin = _fmesh.vertex(0).coord(1);
//double ymax = _fmesh.vertex(0).coord(1);
//for (int i = 1; i < _fmesh.n_vertices(); ++i)
//{
//  const double x = _fmesh.vertex(i).coord(0);
//  const double y = _fmesh.vertex(i).coord(1);
//  if (x < xmin) xmin = x;
//  if (x > xmax) xmax = x;
//  if (y < ymin) ymin = y;
//  if (y > ymax) ymax = y;
//}

//// find nodes that lie on these limits
//std::vector<int> b_nodes;
//for (int i = 0; i < _fmesh.n_vertices(); ++i)
//{
//  const double x = _fmesh.vertex(i).coord(0);
//  const double y = _fmesh.vertex(i).coord(1);
//  if ((fabs(x - xmin) < 1e-14 || fabs(x - xmax) < 1e-14 ||
//       fabs(y - ymin) < 1e-14 || fabs(y - ymax) < 1e-14) &&
//      find(b_nodes.begin(), b_nodes.end(), i) == b_nodes.end())
//    b_nodes.push_back(i);
//}

//#if defined(DEBUG)
//std::cout << "boundary nodes" << std::endl;
//for (int i = 0; i < b_nodes.size(); ++i)
//  std::cout << b_nodes[i] << " ";
//std::cout << std::endl;
//#endif

//// impose Dirichlet boundary condition
//// with ones on diagonal
//MatZeroRows(dense_stiff_mat, b_nodes.size(), &b_nodes[0], 1., solution, system_rhs); // change the matrix
////MatZeroRowsColumns(dense_stiff_mat, b_nodes.size(), &b_nodes[0], 1., solution, system_rhs);
//for (int i = 0; i < b_nodes.size(); ++i)
//  VecSetValue(system_rhs, b_nodes[i], boundary_function(_fmesh.vertex(b_nodes[i])), INSERT_VALUES); // change the rhs vector


//// solve the SLAE
//KSP ksp;
//KSPCreate(PETSC_COMM_WORLD, &ksp);
//KSPSetOperators(ksp, dense_stiff_mat, dense_stiff_mat, SAME_PRECONDITIONER);
//KSPSolve(ksp, system_rhs, solution);

//// check solution
//Vec exact_solution;
//VecDuplicate(solution, &exact_solution);
//for (int i = 0; i < _fmesh.n_vertices(); ++i)
//{
//  Point vert = _fmesh.vertex(i);
//  VecSetValue(exact_solution, i, analytic_solution(vert), INSERT_VALUES);
//}

////  std::cout << "error = " << rel_error(solution, exact_solution) << std::endl;

//Result res(&dof_handler);
//std::string fname = "results.vtu";
//res.write_vtu(fname, solution, exact_solution);



////  double sum_diff = 0.;
////  double sum_an = 0.;
////  for (int i = 0; i < _fmesh.n_vertices(); ++i)
////  {
////    Point vert = _fmesh.vertex(i);
////    std::cout << analytic_solution(vert) << " " << solution[i] << std::endl;
////    sum_diff += (analytic_solution(vert) - solution[i]) * (analytic_solution(vert) - solution[i]);
////    sum_an += analytic_solution(vert) * analytic_solution(vert);
////  }
////  std::cout << "\n\nrelative error = " << sqrt(sum_diff / sum_an) << std::endl;


////  int* nnz = new int[csr_pattern.order()];
////  for (int i = 0; i < csr_pattern.order(); ++i)
////    nnz[i] = csr_pattern.row(i + 1) - csr_pattern.row(i);
////  const int *nnz = csr_pattern.nnz();
////  for (int i = 0; i < csr_pattern.order(); ++i)
////    std::cout << nnz[i] << std::endl;

////  Mat sparse_mat_stiff, sparse_mat_mass;
////  MatCreateSeqAIJ(PETSC_COMM_SELF, csr_pattern.order(), csr_pattern.order(), 0, csr_pattern.nnz(), &sparse_mat_mass);
////  MatCreateSeqAIJ(PETSC_COMM_SELF, csr_pattern.order(), csr_pattern.order(), 0, csr_pattern.nnz(), &sparse_mat_stiff);

////  VecSet(system_rhs, 0.);

////  for (int cell = 0; cell < _fmesh.n_triangles(); ++cell)
////  {
////    Triangle triangle = _fmesh.triangle(cell);
////    triangle.local_mass_matrix(local_mass_mat);
////    triangle.local_stiffness_matrix(local_stiff_mat);
////    triangle.local_rhs_vector(local_rhs_vec, rhs_function, _fmesh.vertices());

////    for (int i = 0; i < triangle.n_dofs(); ++i)
////    {
////      const int dof_i = triangle.dof(i);
////      VecSetValue(system_rhs, dof_i, local_rhs_vec[i], ADD_VALUES);
////      for (int j = 0; j < triangle.n_dofs(); ++j)
////      {
////        const int dof_j = triangle.dof(j);
////        MatSetValue(sparse_mat_mass, dof_i, dof_j, local_mass_mat[i][j], ADD_VALUES);
////        MatSetValue(sparse_mat_stiff, dof_i, dof_j, local_stiff_mat[i][j], ADD_VALUES);
////      }
////    }
////  }

////  MatAssemblyBegin(sparse_mat_mass, MAT_FINAL_ASSEMBLY);
////  MatAssemblyEnd(sparse_mat_mass, MAT_FINAL_ASSEMBLY);

////  MatAssemblyBegin(sparse_mat_stiff, MAT_FINAL_ASSEMBLY);
////  MatAssemblyEnd(sparse_mat_stiff, MAT_FINAL_ASSEMBLY);


////MatView(sparse_mat, PETSC_VIEWER_STDOUT_SELF);
////MatDestroy(&sparse_mat);

////  VecView(system_rhs, PETSC_VIEWER_STDOUT_SELF);
////  VecView(solution, PETSC_VIEWER_STDOUT_SELF);
////  std::cout << "analytic\n";
////  for (int i = 0; i < _fmesh.n_vertices(); ++i)
////  {
////    Point vert = _fmesh.vertex(i);
////    std::cout << analytic_solution(vert) << std::endl;
////  }

////  MatView(dense_stiff_mat, PETSC_VIEWER_STDOUT_SELF);

////  MatDestroy(&dense_mass_mat);
////  MatDestroy(&dense_stiff_mat);

////  VecDestroy(&system_rhs);
////  VecDestroy(&solution);
////  VecDestroy(&exact_solution);



void Acoustic2D::find_bound_nodes(std::vector<int> &b_nodes) const
{
  // find max and min values to know the limits of computational domain.
  // this approach works only if you have a rectangular domain
  double xmin = _fmesh.vertex(0).coord(0);
  double xmax = _fmesh.vertex(0).coord(0);
  double ymin = _fmesh.vertex(0).coord(1);
  double ymax = _fmesh.vertex(0).coord(1);
  for (int i = 1; i < _fmesh.n_vertices(); ++i)
  {
    const double x = _fmesh.vertex(i).coord(0);
    const double y = _fmesh.vertex(i).coord(1);
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
  }

  // find nodes that lie on these limits
  const double tol = 1e-14;
  for (int i = 0; i < _fmesh.n_vertices(); ++i)
  {
    const double x = _fmesh.vertex(i).coord(0);
    const double y = _fmesh.vertex(i).coord(1);
    if ((fabs(x - xmin) < tol || fabs(x - xmax) < tol ||
         fabs(y - ymin) < tol || fabs(y - ymax) < tol) &&
        find(b_nodes.begin(), b_nodes.end(), i) == b_nodes.end())
      b_nodes.push_back(i);
  }
}
