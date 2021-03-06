#include "acoustic2d.h"
#include "parameters.h"
#include "petscksp.h"
#include "fem/auxiliary_functions.h"
#include "fem/finite_element.h"
#include "analytic_functions.h"
#include "fem/result.h"
#include "fem/math_functions.h"
#include "layer.h"
#include "block_of_layers.h"
#include <iostream>
#include <algorithm>
#include <fstream>

using namespace fem;

Acoustic2D::Acoustic2D(Parameters *param)
  : _param(param)
{
  require(_param->FE_ORDER == 1, "This fe order hasn't been implemented");
  require(_param->RES_DIR != "", "Computation environment was not established through Parameter function");
}


Acoustic2D::~Acoustic2D()
{
  //VecDestroy(&_global_rhs);
  //MatDestroy(&_global_mass_mat);
  //MatDestroy(&_global_stiff_mat);
}



void Acoustic2D::solve_rectangles()
{
  require(_param->FE_ORDER == 1, "This fe order is not implemented (" + d2s(_param->FE_ORDER) + ")");

  // create rectangular grid according to parameters defined in _param
  _fmesh.create_rectangular_grid(_param->X_BEG, _param->X_END,
                                 _param->Y_BEG, _param->Y_END,
                                 _param->N_FINE_X, _param->N_FINE_Y);

#if defined(DEBUG)
  std::cout << "n_vertices = " << _fmesh.n_vertices() << std::endl;
  std::cout << "n_rectangles = " << _fmesh.n_rectangles() << std::endl;
#endif

  FiniteElement fe(_param->FE_ORDER);

  DoFHandler dof_handler(&_fmesh);
  dof_handler.distribute_dofs(fe, CG);
#if defined(DEBUG)
  std::cout << "n_dofs = " << dof_handler.n_dofs() << std::endl;
#endif

  // create sparse format based on the distribution of degrees of freedom.
  // since we use first order basis functions, and then
  // all dofs are associated with the mesh vertices,
  // sparse format is based on connectivity of the mesh vertices
  CSRPattern csr_pattern;
  csr_pattern.make_sparse_format(dof_handler, CG);
#if defined(DEBUG)
  std::cout << "csr_order = " << csr_pattern.order() << std::endl;
#endif

  expect(csr_pattern.order() == dof_handler.n_dofs(), "Error");


  // allocate memory
  VecCreateSeq(PETSC_COMM_SELF, csr_pattern.order(), &_global_rhs);
  MatCreateSeqAIJ(PETSC_COMM_WORLD, csr_pattern.order(), csr_pattern.order(), 0, csr_pattern.nnz(), &_global_mass_mat);
  MatCreateSeqAIJ(PETSC_COMM_WORLD, csr_pattern.order(), csr_pattern.order(), 0, csr_pattern.nnz(), &_global_stiff_mat);

  // allocate the memory for local matrices and vectors
  double **local_mass_mat = new double*[Rectangle::n_dofs_first];
  double **local_stiff_mat = new double*[Rectangle::n_dofs_first];
  for (unsigned int i = 0; i < Rectangle::n_dofs_first; ++i)
  {
    local_mass_mat[i] = new double[Rectangle::n_dofs_first];
    local_stiff_mat[i] = new double[Rectangle::n_dofs_first];
  }

  // fill up the array of coefficients alpha and beta
  if (_param->CREATE_BIN_LAYERS_FILE)
  {
    if (_param->WHAT_BIN_LAYERS_FILE == "3")
      create_3_bin_layers_file();
    else if (_param->WHAT_BIN_LAYERS_FILE == "slop")
      create_slop_bin_layers_file();
    else
      require(false, "Unknown sort of bin layer file to be created: " + _param->WHAT_BIN_LAYERS_FILE);
  }
//  if (_param->CREATE_AVE_LAYERS_FILE)
//    create_ave_layers_file();

  if (_param->COEF_SAVED_PER_VERT) // if coefficients have been saved in a file,
  {
    import_coefficients_per_vertex(_param->COEF_FILE); // we just read them and convert from vertex-wise distribution to cell-wise one
  }
  else // if the coefficients haven't been saved before, we distribute them according to layers file or other ways
  {
    if (_param->USE_LAYERS_FILE)
      coefficients_initialization();
    else
    {
      _coef_alpha.resize(_fmesh.n_rectangles(), _param->COEF_A_VALUES[0]);
      _coef_beta.resize(_fmesh.n_rectangles(), _param->COEF_B_VALUES[0]);
      for (unsigned int cell = 0; cell < _fmesh.n_rectangles(); ++cell)
      {
        if ((int)_fmesh.rectangle(cell).material_id() == _param->INCL_DOMAIN)
        {
          _coef_alpha[cell] = _param->COEF_A_VALUES[1]; // coefficient in the inclusion
          _coef_beta[cell]  = _param->COEF_B_VALUES[1];  // coefficient in the inclusion
        }
      }
    }
    require(!(_param->SAVE_COEF_PER_CELL && _param->SAVE_COEF_PER_VERT), "Conflicting options");
    if (_param->SAVE_COEF_PER_CELL) // if we need to save coefficients after distribution - one coef per cell,
      export_coefficients_per_cell(_param->COEF_FILE); // we just save them into a file
    else if (_param->SAVE_COEF_PER_VERT) // if we need to save coefficients after distribution - one per vertex,
      export_coefficients_per_vertex(_param->COEF_FILE); // we convert them to vertex-wise format and export them into a file
  }



  // assemble the matrices
  for (unsigned int cell = 0; cell < _fmesh.n_rectangles(); ++cell)
  {
    const Rectangle rectangle = _fmesh.rectangle(cell);
    rectangle.local_mass_matrix(_coef_alpha[cell], local_mass_mat);
    rectangle.local_stiffness_matrix(_coef_beta[cell], local_stiff_mat);

    for (unsigned int i = 0; i < rectangle.n_dofs(); ++i)
    {
      const unsigned int dof_i = rectangle.dof(i);
      for (unsigned int j = 0; j < rectangle.n_dofs(); ++j)
      {
        const unsigned int dof_j = rectangle.dof(j);
        MatSetValue(_global_mass_mat, dof_i, dof_j, local_mass_mat[i][j], ADD_VALUES);
        MatSetValue(_global_stiff_mat, dof_i, dof_j, local_stiff_mat[i][j], ADD_VALUES);
      }
    }
  }

  // free the memory
  for (unsigned int i = 0; i < Rectangle::n_dofs_first; ++i)
  {
    delete[] local_mass_mat[i];
    delete[] local_stiff_mat[i];
  }
  delete[] local_mass_mat;
  delete[] local_stiff_mat;

  MatAssemblyBegin(_global_mass_mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(_global_mass_mat, MAT_FINAL_ASSEMBLY);

  MatAssemblyBegin(_global_stiff_mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(_global_stiff_mat, MAT_FINAL_ASSEMBLY);

  if (_param->TIME_SCHEME == EXPLICIT)
    solve_explicit_rectangles(dof_handler, csr_pattern);
  else if(_param->TIME_SCHEME == CRANK_NICOLSON)
    solve_crank_nicolson(dof_handler, csr_pattern);
  else
    require(false, "Unknown time discretization scheme");
}



void Acoustic2D::solve_explicit_rectangles(const DoFHandler &dof_handler, const CSRPattern &csr_pattern)
{
  require(_param->FE_ORDER == 1, "This fe order is not implemented (" + d2s(_param->FE_ORDER) + ")");

  // create vectors
  Vec solution; // numerical solution on the current (n-th) time step
  Vec solution_1; // numerical solution on the previous ((n-1)-th) time step
  Vec solution_2; // numerical solution on the preprevious ((n-2)-th) time step
  //Vec exact_solution; // analytic solution

  VecDuplicate(_global_rhs, &solution);
  VecDuplicate(_global_rhs, &solution_1);
  VecDuplicate(_global_rhs, &solution_2);
  //VecDuplicate(_global_rhs, &exact_solution);

  // for brevity
  const double dt = _param->TIME_STEP;

  // fill vectors with solution on the 0-th and 1-st time steps
  const InitialSolution init_solution;
  for (unsigned int d = 0; d < dof_handler.n_dofs(); ++d)
  {
    VecSetValue(solution_2, d, init_solution.value(dof_handler.dof(d), _param->TIME_BEG), INSERT_VALUES);
    VecSetValue(solution_1, d, init_solution.value(dof_handler.dof(d), _param->TIME_BEG + dt), INSERT_VALUES);
  }

  // make a SLAE rhs vector
  Vec system_rhs, temp;
  VecCreateSeq(PETSC_COMM_SELF, csr_pattern.order(), &system_rhs);
  VecDuplicate(system_rhs, &temp);

  // boundary nodes
  const std::vector<int> &b_nodes = _fmesh.boundary_vertices();

  // system matrix equal to global mass matrix
  Mat system_mat;
  MatConvert(_global_mass_mat, MATSAME, MAT_INITIAL_MATRIX, &system_mat); // allocate memory and copy values

  // impose Dirichlet boundary condition
  // with ones on diagonal
  MatZeroRows(system_mat, b_nodes.size(), &b_nodes[0], 1., solution, system_rhs); // change the matrix

  // SLAE solver
  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, system_mat, system_mat, SAME_PRECONDITIONER);
  KSPSetTolerances(ksp, 1e-12, 1e-30, 1e+5, 10000);

  double *local_rhs_vec = new double[Rectangle::n_dofs_first];

  require(_param->N_TIME_STEPS > 1, "There is no time steps to perform: n_time_steps = " + d2s(_param->N_TIME_STEPS));

  const RHSFunction rhs_function(*_param);

  if (_param->PRINT_INFO)
    std::cout << "time loop started..." << std::endl;

  for (unsigned int time_step = 2; time_step <= _param->N_TIME_STEPS; ++time_step)
  {
    const double time = _param->TIME_BEG + time_step * dt; // current time
    VecSet(system_rhs, 0.); // zeroing the system rhs vector
    VecSet(solution, 0.); // zeroing the solution vector

    // assemble some parts of system rhs vector
    for (unsigned int cell = 0; cell < _fmesh.n_rectangles(); ++cell)
    {
      Rectangle rectangle = _fmesh.rectangle(cell);
      rectangle.local_rhs_vector(rhs_function, dof_handler.dofs(), time - dt, local_rhs_vec); // rhs function on the previous time step
#if defined(DEBUG)
      std::cout << "\trectangle " << cell << "\n\t";
      for (int i = 0; i < rectangle.n_dofs(); ++i)
        std::cout << local_rhs_vec[i] << " ";
      std::cout << std::endl;
#endif
      for (unsigned int i = 0; i < rectangle.n_dofs(); ++i)
      {
        const unsigned int dof_i = rectangle.dof(i);
        VecSetValue(system_rhs, dof_i, local_rhs_vec[i], ADD_VALUES);
      }
    } // rhs part assembling

    MatMult(_global_stiff_mat, solution_1, temp);
    VecAXPBY(system_rhs, -dt*dt, dt*dt, temp);

    MatMult(_global_mass_mat, solution_2, temp);
    VecAXPY(system_rhs, -1., temp);

    MatMult(_global_mass_mat, solution_1, temp);
    VecAXPY(system_rhs, 2., temp);

    // impose Dirichlet boundary condition
    const BoundaryFunction boundary_function;
    for (unsigned int i = 0; i < b_nodes.size(); ++i)
      VecSetValue(system_rhs, b_nodes[i], boundary_function.value(_fmesh.vertex(b_nodes[i]), time), INSERT_VALUES); // change the rhs vector

    // solve the SLAE
    KSPSolve(ksp, system_rhs, solution);

    // reassign the solutions on the previuos time steps
    VecCopy(solution_1, solution_2);
    VecCopy(solution,   solution_1);

    // check solution
//    for (int i = 0; i < _fmesh.n_vertices(); ++i)
//    {
//      Point vert = _fmesh.vertex(i);
//      VecSetValue(exact_solution, i, an_solution(vert, time), INSERT_VALUES);
//    }
//    std::cout << "time step = " << time_step << " time = " << time << " relative error = " << rel_error(solution, exact_solution) << std::endl;

    if ((_param->PRINT_VTU && (time_step % _param->VTU_STEP == 0)) || (time_step == _param->N_TIME_STEPS))
    {
      Result res(&dof_handler);
      std::string fname = _param->VTU_DIR + "/res-" + d2s(time_step) + ".vts";
      if (time_step == _param->N_TIME_STEPS && _param->EXPORT_COEFFICIENTS) // we export coefficients only for the last step
        res.write_vts(fname, _param->N_FINE_X, _param->N_FINE_Y, solution, NULL, _coef_alpha, _coef_beta);
      else
        res.write_vts(fname, _param->N_FINE_X, _param->N_FINE_Y, solution); //, exact_solution, _coef_alpha, _coef_beta);
    }

    if (_param->PRINT_INFO)
    {
      double rhs_norm;
      VecNorm(system_rhs, NORM_2, &rhs_norm);
      double system_rhs_norm;
      VecNorm(system_rhs, NORM_2, &system_rhs_norm);
      double norm;
      VecNorm(solution, NORM_2, &norm);
      std::cout.setf(std::ios::scientific);
      std::cout.precision(4);
      std::cout << "  step " << time_step << " norm " << norm << " rhs_norm " << rhs_norm
                << " sys_rhs_norm " << system_rhs_norm << std::endl;
    }

    if ((_param->SAVE_SOL && (time_step % _param->SOL_STEP == 0)) || (time_step == _param->N_TIME_STEPS))
    {
      // extract data from PETSc vector
      std::vector<int> idx(csr_pattern.order());
      std::iota(idx.begin(), idx.end(), 0); // idx = { 0, 1, 2, 3, .... }
      std::vector<double> solution_values(csr_pattern.order());
      VecGetValues(solution, csr_pattern.order(), &idx[0], &solution_values[0]);
      // write the solution to the file
      std::string fname = _param->SOL_DIR + "/sol-" + d2s(time_step) + ".dat";
      std::ofstream out(fname);
      out.setf(std::ios::scientific);
      out.precision(16);
      //out << csr_pattern.order() << "\n";
      for (unsigned int i = 0; i < csr_pattern.order(); ++i)
        out << solution_values[i] << "\n";
      out.close();
    }

  } // time loop

  KSPDestroy(&ksp);

  MatDestroy(&system_mat);

  delete[] local_rhs_vec;

  VecDestroy(&solution);
  VecDestroy(&solution_1);
  VecDestroy(&solution_2);
  //VecDestroy(&exact_solution);
  VecDestroy(&system_rhs);
  VecDestroy(&temp);
}



void Acoustic2D::coefficients_initialization()
{
  std::ifstream in(_param->LAYERS_FILE.c_str());
  require(in, "File " + _param->LAYERS_FILE + " cannot be opened");

  const std::vector<Rectangle> &cells = _fmesh.rectangles(); // all mesh cells
  _coef_alpha.resize(cells.size(), 0); // coefficient alpha in each cell is 0 by default
  _coef_beta.resize(cells.size(),  0); // coefficient beta in each cell is 0 by default

  // since the layers are distributed horisontally (or nearly horisontally) in most cases
  // the thickness of the layer is associated with the vertical axis
  // and expressed in percents according to the height (depth) of the domain
  const double Hy = _fmesh.max_coord().coord(1) - _fmesh.min_coord().coord(1);

  // general structure of the layers file is the following:
  // M                        (M - number of blocks with different layers distribution (different angle, or amount, for example))
  // v1 v2 N angle            (v1 - relative heigth of the beginning of the layer,
  //                           v2 - realtive height of its ending,
  //                           N - the number of layers inside one part of the domain,
  //                           angle - slope angle of all layers inside this block)
  // h1 a1 b1                 (h1 relative heigth of the layer 1 within this block, a1, b1 - coefficients)
  // h2 a2 b2                 (h2 relative heigth of the layer 2 within this block, a2, b2 - coefficients)
  // ...
  // hN aN bN                 (hN relative heigth of the layer N within this block, aN, bN - coefficients)
  // v1 v2 N angle            ANOTHER block
  // and so on

  unsigned int n_blocks; // the number of block of layers
  in >> n_blocks;
  require(n_blocks > 0, "The number of blocks is 0");
  std::vector<BlockOfLayers> blocks(n_blocks); // allocate the memory for all blocks
  std::vector<double> block_beg(n_blocks), block_end(n_blocks); // in percents

  for (unsigned int bl = 0; bl < n_blocks; ++bl)
  {
    double layer_angle; // slope angle of the layers
    unsigned int n_layers; // the number of layers
    in >> block_beg[bl] >> block_end[bl] >> n_layers >> layer_angle;
    require(block_beg[bl] >= 0 &&
            block_end[bl] <= 100 &&
            block_beg[bl] < block_end[bl], "The height of layer is wrong");

    const double right_angle = 90.; // right angle in degrees
    require(fabs(layer_angle) < right_angle, "Angle doesn't belong to correct range: (-90, 90), its value : " + d2s(layer_angle));
    require(fabs(fabs(layer_angle) - right_angle) > math::FLOAT_NUMBERS_EQUALITY_TOLERANCE, "Angle is equal to right angle (90), what is prohibited");

    require(n_layers > 0, "The number of layers is 0");
    std::vector<double> layer_h_percent(n_layers); // the thicknesses of the layers in percents
    std::vector<double> layer_coef_alpha(n_layers); // the coefficients alpha in each layer
    std::vector<double> layer_coef_beta(n_layers); // the coefficients beta in each layer

    const double x0 = _fmesh.min_coord().coord(0);
    const double x1 = _fmesh.max_coord().coord(0);
    double y0, y1;
    if (fabs(block_beg[bl]) < math::FLOAT_NUMBERS_EQUALITY_TOLERANCE) // block_beg is 0
      y0 = _fmesh.min_coord().coord(1);
    else
      y0 = _fmesh.min_coord().coord(1) + 0.01 * block_beg[bl] * Hy;
    if (fabs(block_end[bl] - 100.) < math::FLOAT_NUMBERS_EQUALITY_TOLERANCE) // block_end is 100
      y1 = _fmesh.max_coord().coord(1);
    else
      y1 = _fmesh.min_coord().coord(1) + 0.01 * block_end[bl] * Hy;

    // block limits
    const Point min_block_point(x0, y0);
    const Point max_block_point(x1, y1);

    // read the info abount thickness and coefficients of each layer
    double total_h_percent = 0.; // in the end of the day total_h should be 100 (because it's in percent). in other case those thicknesses don't make sense
    for (unsigned int i = 0; i < n_layers; ++i)
    {
      // read info
      in >> layer_h_percent[i]
         >> layer_coef_alpha[i]
         >> layer_coef_beta[i];
      total_h_percent += layer_h_percent[i];
    }
    require(fabs(total_h_percent - 100.) < math::FLOAT_NUMBERS_EQUALITY_REDUCED_TOLERANCE,
            "The total thickness of all layers is not 100%. It is " + d2s(total_h_percent, true, 14));

    // block initialization
    blocks[bl].init(min_block_point, max_block_point,
                    n_layers, layer_angle,
                    layer_h_percent,
                    layer_coef_alpha,
                    layer_coef_beta);

  } // for all blocks

  // check that blocks don't intersect each other and fill up the whole domain
  double total_block_h = 0.; // in percent
  for (unsigned int i = 0; i < n_blocks; ++i)
    total_block_h += block_end[i] - block_beg[i];
  require(fabs(total_block_h - 100.) < math::FLOAT_NUMBERS_EQUALITY_TOLERANCE,
          "The blocks either intersect each other or don't fill up the whole domain");

  // distribute the coefficients in each cell according to the layers
  for (unsigned int i = 0; i < cells.size(); ++i)
  {
    bool coef_found = false;
    for (unsigned int j = 0; j < n_blocks && coef_found == false; ++j)
    {
      if (blocks[j].contains_element(cells[i], _fmesh.vertices()))
      {
        blocks[j].get_coefs(cells[i], _fmesh.vertices(),
                            _coef_alpha[i], _coef_beta[i]);
        coef_found = true;
      }
    }
    require(coef_found, "The cell number " + d2s(i) + " doesn't belong to any block");
  }

  if (_param->USE_AVERAGED) // if we use averaged coefficient on a part of a domain
  {
    const double Hy = _fmesh.max_coord().coord(1) - _fmesh.min_coord().coord(1); // the height (depth) of the domain
    const double thickness_limit = 10;

    // we average the coefficients alpha and beta on those parts that occupy less than 10% of the domain
    for (unsigned int j = 0; j < n_blocks; ++j)
    {
      // averaged coefficients can be different for each block
      double aver_alpha = 0;
      double aver_beta  = 0;
      double total_mes  = 0; // the sum of the thicknesses of all layers which coefficients we average

      if (blocks[j].n_layers() > 1)
      {
        for (unsigned int i = 0; i < cells.size(); ++i)
        {
          if (blocks[j].contains_element(cells[i], _fmesh.vertices()))
          {
            const Layer &layer = blocks[j].layer_which_contains(cells[i], _fmesh.vertices()); // find the layer which contains the cell
            const double relative_thickness = layer.thickness() / Hy * 100; // the thickness of the layer in respect with the height of the whole domain in percent
            if (relative_thickness < thickness_limit) // if the relative thickness is less than this number of percent we average the coefficient, otherwise we don't
            {
              const double cell_mes = cells[i].mes(); // the measure (area, volume) of the cell
              aver_alpha += cell_mes / _coef_alpha[i];
              aver_beta  += cell_mes / _coef_beta[i];
              total_mes  += cell_mes;
              _coef_alpha[i] = -1; // we do that to know the coefficients in which cells will have to be renewed
              _coef_beta[i]  = -1;
            }
          }
        } // for each cell

        aver_alpha /= total_mes;
        aver_beta  /= total_mes;

        aver_alpha = 1. / aver_alpha;
        aver_beta  = 1. / aver_beta;

        // now we distribute averaged coefficients for all cells which were changed
        for (unsigned int i = 0; i < cells.size(); ++i)
        {
          if (fabs(_coef_alpha[i] + 1) < math::FLOAT_NUMBERS_EQUALITY_TOLERANCE)
          {
            _coef_alpha[i] = aver_alpha;
            _coef_beta[i]  = aver_beta;
          }
        }

      } // if a block has more than 1 layer
    } // for each block
  } // if we use averaged coefficients
}



void Acoustic2D::create_3_bin_layers_file() const
{
  std::string fname = _param->LAYERS_DIR + "/lay_3_bin_" + d2s(_param->H_BIN_LAYER_PERCENT) + "_" + _param->LAYERS_FILE_SUFFIX + ".dat";

  std::cout << "\ncreation of binary layers file: " << fname << std::endl;
  std::ofstream out(fname.c_str());
  require(out, "Cannot open file " + fname);

  out.setf(std::ios::scientific);
  out.precision(14);

  const unsigned int n_blocks = 3;
  out << n_blocks << "\n";

  const unsigned int block_beg[] = { 0, 80, 20 };
  const unsigned int block_end[] = { 20, 100, 80 };
  //const unsigned int n_layers[]  = { 1, 1, 20 };

  if (fabs(block_end[0] - block_beg[0]) > math::FLOAT_NUMBERS_EQUALITY_TOLERANCE)
  {
    out << block_beg[0] << " " << block_end[0] << " 1 0\n";
    out << "100 " << _param->COEF_A_VALUES[0] << " " << _param->COEF_B_VALUES[0] << "\n";
  }

  if (fabs(block_end[1] - block_beg[1]) > math::FLOAT_NUMBERS_EQUALITY_TOLERANCE)
  {
    out << block_beg[1] << " " << block_end[1] << " 1 0\n";
    out << "100 " << _param->COEF_A_VALUES[0] << " " << _param->COEF_B_VALUES[0] << "\n";
  }

  const double Hy = _fmesh.max_coord().coord(1) - _fmesh.min_coord().coord(1);
  const double thickness = 0.01 * (block_end[2] - block_beg[2]) * Hy;
  const double h_thin_layer = 0.01 * _param->H_BIN_LAYER_PERCENT * Hy; // some meters if the domain is 1000 x 1000 m
  const unsigned int n_thin_layers = (unsigned)(thickness / h_thin_layer);
  double coef_a, coef_b; // for brevity

  out << block_beg[2] << " " << block_end[2] << " " << n_thin_layers << " 0\n";
  for (unsigned int i = 0; i < n_thin_layers; ++i)
  {
    // for even layers we use _param->COEF_*_VALUES[1]
    // for odd  layers we use _param->COEF_*_VALUES[0]
    if (i % 2 == 0)
    {
      coef_a = _param->COEF_A_VALUES[1];
      coef_b = _param->COEF_B_VALUES[1];
    }
    else
    {
      coef_a = _param->COEF_A_VALUES[0];
      coef_b = _param->COEF_B_VALUES[0];
    }

    out << 100. / (double)n_thin_layers << " "
        << coef_a << " "
        << coef_b << "\n"; // binary medium
  }

#if defined(DEBUG)
  std::cout << "thickness = " << thickness << std::endl;
  std::cout << "h_thin_layer = " << h_thin_layer << std::endl;
  std::cout << "n_thin_layers = " << n_thin_layers << std::endl;
#endif

  out.close();
}



void Acoustic2D::create_slop_bin_layers_file() const
{
  std::string fname = _param->LAYERS_DIR + "/lay_slop_bin_" + d2s(_param->H_BIN_LAYER_PERCENT) + "_" + _param->LAYERS_FILE_SUFFIX + ".dat";

  std::cout << "\ncreation of binary layers file: " << fname << std::endl;
  std::ofstream out(fname.c_str());
  require(out, "Cannot open file " + fname);

  out.setf(std::ios::scientific);
  out.precision(14);

  const unsigned int n_blocks = 1;
  out << n_blocks << "\n";

  const unsigned int block_beg[] = { 0, 0, 0 }; // { 0, 80, 20 };
  const unsigned int block_end[] = { 0, 0, 100 }; //{ 20, 100, 80 };
  //const unsigned int n_layers[]  = { 1, 1, 20 };

  if (fabs(block_end[0] - block_beg[0]) > math::FLOAT_NUMBERS_EQUALITY_TOLERANCE)
  {
    out << block_beg[0] << " " << block_end[0] << " 1 0\n";
    out << "100 " << _param->COEF_A_VALUES[0] << " " << _param->COEF_B_VALUES[0] << "\n";
  }

  if (fabs(block_end[1] - block_beg[1]) > math::FLOAT_NUMBERS_EQUALITY_TOLERANCE)
  {
    out << block_beg[1] << " " << block_end[1] << " 1 0\n";
    out << "100 " << _param->COEF_A_VALUES[0] << " " << _param->COEF_B_VALUES[0] << "\n";
  }

  const double Hy = _fmesh.max_coord().coord(1) - _fmesh.min_coord().coord(1);
  const double thickness = 0.01 * (block_end[2] - block_beg[2]) * Hy;
  const double h_thin_layer = 0.01 * _param->H_BIN_LAYER_PERCENT * Hy; // some meters if the domain is 1000 x 1000 m
  const double thick_layer = 20; // we have 2 thick sloped layers
  const unsigned int n_thin_layers = (unsigned)(thickness / h_thin_layer) - 2;
  double coef_a, coef_b; // for brevity

  out << block_beg[2] << " " << block_end[2] << " " << n_thin_layers + 2 << " 45\n";
  out << thick_layer << " " << _param->COEF_A_VALUES[0] << " " << _param->COEF_B_VALUES[0] << "\n";
  for (unsigned int i = 0; i < n_thin_layers; ++i)
  {
    // for even layers we use _param->COEF_*_VALUES[1]
    // for odd  layers we use _param->COEF_*_VALUES[0]
    if (i % 2 == 0)
    {
      coef_a = _param->COEF_A_VALUES[1];
      coef_b = _param->COEF_B_VALUES[1];
    }
    else
    {
      coef_a = _param->COEF_A_VALUES[0];
      coef_b = _param->COEF_B_VALUES[0];
    }

    out << (100. - 2 * thick_layer) / (double)n_thin_layers << " "
        << coef_a << " "
        << coef_b << "\n"; // binary medium
  }
  out << thick_layer << " " << _param->COEF_A_VALUES[0] << " " << _param->COEF_B_VALUES[0] << "\n";

#if defined(DEBUG)
  std::cout << "thickness = " << thickness << std::endl;
  std::cout << "h_thin_layer = " << h_thin_layer << std::endl;
  std::cout << "n_thin_layers = " << n_thin_layers << std::endl;
#endif

  out.close();
}



//void Acoustic2D::create_ave_layers_file() const
//{
//  std::string fname = _param->LAYERS_DIR + "/lay_3_ave_" + d2s(_param->H_BIN_LAYER_PERCENT) + "_" + _param->LAYERS_FILE_SUFFIX + ".dat";

//  std::cout << "\ncreation of average layers file: " << fname << std::endl;
//  std::ofstream out(fname.c_str());
//  require(out, "Cannot open file " + fname);

//  out.setf(std::ios::scientific);
//  out.precision(14);

//  const unsigned int n_blocks = 3;
//  out << n_blocks << "\n";

//  const unsigned int block_beg[] = { 0, 80, 20 };
//  const unsigned int block_end[] = { 20, 100, 80 };
//  //const unsigned int n_layers[]  = { 1, 1, 20 };

//  if (fabs(block_end[0] - block_beg[0]) > math::FLOAT_NUMBERS_EQUALITY_TOLERANCE)
//  {
//    out << block_beg[0] << " " << block_end[0] << " 1 0\n";
//    out << "100 " << _param->COEF_A_VALUES[0] << " " << _param->COEF_B_VALUES[0] << "\n";
//  }

//  if (fabs(block_end[1] - block_beg[1]) > math::FLOAT_NUMBERS_EQUALITY_TOLERANCE)
//  {
//    out << block_beg[1] << " " << block_end[1] << " 1 0\n";
//    out << "100 " << _param->COEF_A_VALUES[0] << " " << _param->COEF_B_VALUES[0] << "\n";
//  }

//  const double Hy = _fmesh.max_coord().coord(1) - _fmesh.min_coord().coord(1);
//  const double thickness = 0.01 * (block_end[2] - block_beg[2]) * Hy;
//  const double h_thin_layer = 0.01 * _param->H_BIN_LAYER_PERCENT * Hy; // some meters if the domain is 1000 x 1000 m
//  const unsigned int n_thin_layers = (unsigned)(thickness / h_thin_layer);

//  double coef_a, coef_b; // for brevity
//  double coef_a_aver = 0.;
//  double coef_b_aver = 0.;

//  for (unsigned int i = 0; i < n_thin_layers; ++i)
//  {
//    // for even layers we use _param->COEF_*_VALUES[1]
//    // for odd  layers we use _param->COEF_*_VALUES[0]
//    if (i % 2 == 0)
//    {
//      coef_a = _param->COEF_A_VALUES[1];
//      coef_b = sqrt(_param->COEF_B_VALUES[1]); // we need a value of velocity
//    }
//    else
//    {
//      coef_a = _param->COEF_A_VALUES[0];
//      coef_b = sqrt(_param->COEF_B_VALUES[0]); // we need a value of velocity
//    }

////    coef_a_aver += h_thin_layer * coef_a;
////    coef_b_aver += h_thin_layer * coef_b;
//    coef_a_aver += h_thin_layer / coef_a;
//    coef_b_aver += h_thin_layer / coef_b;
//  }
//  coef_a_aver /= thickness;
//  coef_b_aver /= thickness;

//  coef_a_aver = 1. / coef_a_aver;
//  coef_b_aver = 1. / coef_b_aver; // this is averaged velocity
//  coef_b_aver *= coef_b_aver; // we need velocity squared as a coefficient beta

//#if defined(DEBUG)
//  std::cout << "thickness = " << thickness << std::endl;
//  std::cout << "h_thin_layer = " << h_thin_layer << std::endl;
//  std::cout << "n_thin_layers = " << n_thin_layers << std::endl;
//  std::cout << "coef_a_aver = " << coef_a_aver << std::endl;
//  std::cout << "coef_b_aver = " << coef_b_aver << std::endl;
//#endif

//  out << block_beg[2] << " " << block_end[2] << " 1 0\n";
//  out << "100 " << coef_a_aver << " " << coef_b_aver << "\n"; // average medium

//  out.close();
//}



void Acoustic2D::export_coefficients_per_vertex(const std::string &filename) const
{
  // since a part of software that will use this coefficients distribution
  // is based on suggestion that coefficients are vertex-wise distributed
  // instead of cell-wise distribution accepted in this code,
  // we need to convert our array of coefficients to another one.
  // the rule of conversion is quite simple.
  // assume we have a mesh, and we need to define the value
  // of the suitable coefficient at its vertices.
  // v6 --- v7 --- v8
  // |  c2  |  c3  |
  // v3 --- v4 --- v5
  // |  c0  |  c1  |
  // v0 --- v1 --- v2
  // vi (i=0,..,8) - vertices
  // ci (i=0,..,3) - coefficients
  //
  // we do the following averaging for the vertices of the first cell
  // c(v0) = c0
  // c(v1) = 1/2 * (c0 + c1)
  // c(v3) = 1/2 * (c0 + c2)
  // c(v4) = 1/4 * (c0 + c1 + c2 + c3)
  // that's it.

  const unsigned int nx = _param->N_FINE_X; // for brevity
  const unsigned int ny = _param->N_FINE_Y; // for brevity

  expect(_fmesh.n_rectangles() == nx * ny, "The number of elements is somehow incorrect");
  expect(_fmesh.n_vertices() == (nx + 1) * (ny + 1), "The number of elements is somehow incorrect");

  std::vector<double> coeff_a_vert(_fmesh.n_vertices());
  std::vector<double> coeff_b_vert(_fmesh.n_vertices());

  // left lower corner
  coeff_a_vert[0] = _coef_alpha[0];
  coeff_b_vert[0] = _coef_beta[0];

  // right lower corner
  coeff_a_vert[nx] = _coef_alpha[nx - 1];
  coeff_b_vert[nx] = _coef_beta[nx - 1];

  // left upper corner
  coeff_a_vert[ny * (nx + 1)] = _coef_alpha[(ny - 1) * nx];
  coeff_b_vert[ny * (nx + 1)] = _coef_beta[(ny - 1) * nx];

  // right upper corner
  coeff_a_vert[(ny + 1) * (nx + 1) - 1] = _coef_alpha[ny * nx - 1];
  coeff_b_vert[(ny + 1) * (nx + 1) - 1] = _coef_beta[ny * nx - 1];

  // bottom line except corners (shared by 2 elements)
  for (unsigned int i = 1; i < nx; ++i)
  {
    coeff_a_vert[i] = 0.5 * (_coef_alpha[i - 1] + _coef_alpha[i]);
    coeff_b_vert[i] = 0.5 * (_coef_beta[i - 1] + _coef_beta[i]);
  }

  // top line except corners (shared by 2 elements)
  unsigned int start_v = ny * (nx + 1); // start number for vertices
  unsigned int start_c = (ny - 1) * nx; // start number for cells
  for (unsigned int i = 1; i < nx; ++i)
  {
    coeff_a_vert[start_v + i] = 0.5 * (_coef_alpha[start_c + i - 1] + _coef_alpha[start_c + i]);
    coeff_b_vert[start_v + i] = 0.5 * (_coef_beta[start_c + i - 1] + _coef_beta[start_c + i]);
  }

  // left line except corners (shared by 2 elements)
  //unsigned int start_v = ny * (nx + 1); // start number for vertices
  //unsigned int start_c = (ny - 1) * nx; // start number for cells
  for (unsigned int i = 1; i < ny; ++i)
  {
    coeff_a_vert[i * (nx + 1)] = 0.5 * (_coef_alpha[(i - 1) * nx] + _coef_alpha[i * nx]);
    coeff_b_vert[i * (nx + 1)] = 0.5 * (_coef_beta[(i - 1) * nx] + _coef_beta[i * nx]);
  }

  // right line except corners (shared by 2 elements)
  start_v = nx; // start number for vertices
  start_c = nx - 1; // start number for cells
  for (unsigned int i = 1; i < ny; ++i)
  {
    coeff_a_vert[start_v + i * (nx + 1)] = 0.5 * (_coef_alpha[start_c + (i - 1) * nx] + _coef_alpha[start_c + i * nx]);
    coeff_b_vert[start_v + i * (nx + 1)] = 0.5 * (_coef_beta[start_c + (i - 1) * nx] + _coef_beta[start_c + i * nx]);
  }

  // other vertices (internal - shared by 4 elements)
  for (unsigned int i = 1; i < ny; ++i)
  {
    for (unsigned int j = 1; j < nx; ++j)
    {
      coeff_a_vert[i * (nx + 1) + j] = 0.25 * (_coef_alpha[(i - 1) * nx + j - 1] +
                                               _coef_alpha[(i - 1) * nx + j - 0] +
                                               _coef_alpha[(i - 0) * nx + j - 1] +
                                               _coef_alpha[(i - 0) * nx + j - 0]);
      coeff_b_vert[i * (nx + 1) + j] = 0.25 * (_coef_beta[(i - 1) * nx + j - 1] +
                                               _coef_beta[(i - 1) * nx + j - 0] +
                                               _coef_beta[(i - 0) * nx + j - 1] +
                                               _coef_beta[(i - 0) * nx + j - 0]);
    }
  }

  std::ofstream out(filename.c_str());
  require(out, "File " + filename + " cannot be opened");
  out.setf(std::ios::scientific);
  out.precision(14);

  out << _fmesh.n_vertices() << "\n";
  for (unsigned int i = 0; i < _fmesh.n_vertices(); ++i)
    out << coeff_a_vert[i] << " " << coeff_b_vert[i] << "\n";

  out.close();
}



void Acoustic2D::export_coefficients_per_cell(const std::string &filename) const
{
  // since we already have cell-wise distribution of coefficients
  // we just save them in a file

  const unsigned int nx = _param->N_FINE_X; // for brevity
  const unsigned int ny = _param->N_FINE_Y; // for brevity

  expect(_coef_alpha.size() == nx * ny, "The size of coef alpha vector is somehow incorrect");
  expect(_coef_beta.size()  == nx * ny, "The size of coef beta vector is somehow incorrect");

  std::ofstream out(filename.c_str());
  require(out, "File " + filename + " cannot be opened");
  out.setf(std::ios::scientific);
  out.precision(14);

  out << nx * ny << "\n";
  for (unsigned int i = 0; i < nx * ny; ++i)
    out << _coef_alpha[i] << " " << _coef_beta[i] << "\n";

  out.close();
}



void Acoustic2D::import_coefficients_per_vertex(const std::string &filename)
{
  // for some explanation see comments from export_coefficients function.
  // here we read a set of coefficients and convert them to cell-wise representation,
  // in such a way that if we have a mesh
  // v6 --- v7 --- v8
  // |  c2  |  c3  |
  // v3 --- v4 --- v5
  // |  c0  |  c1  |
  // v0 --- v1 --- v2
  // vi (i=0,..,8) - vertices
  // ci (i=0,..,3) - coefficients
  //
  // then we get the coefficients according to the following rule
  // c0 = 1/4 * (c(v0) + c(v1) + c(v3) + c(v4))
  // c1 = 1/4 * (c(v1) + c(v2) + c(v4) + c(v5))
  // etc.

  std::ifstream in(filename.c_str());
  require(in, "File " + filename + " cannot be opened");

  unsigned int n_values;
  in >> n_values;
  require(n_values == _fmesh.n_vertices(),
          "The number of coefficients from the file " + filename + " doesn't correspond to the current number of mesh vertices");

  std::vector<double> coeff_a_vert(n_values);
  std::vector<double> coeff_b_vert(n_values);

  for (unsigned int i = 0; i < n_values; ++i)
    in >> coeff_a_vert[i] >> coeff_b_vert[i];

  in.close();

  // allocate the memory for cell-wise distributed coefficients
  _coef_alpha.resize(_fmesh.n_rectangles(), 0); // initialized by 0
  _coef_beta.resize(_fmesh.n_rectangles(), 0); // initialized by 0

  // now we distribute the coefficients by cells
  for (unsigned int cell = 0; cell < _fmesh.n_rectangles(); ++cell)
  {
    const Rectangle &element = _fmesh.rectangle(cell);
    for (unsigned int v = 0; v < Rectangle::n_vertices; ++v)
    {
      const unsigned int vert = element.vertex(v); // the number of a current vertex
      _coef_alpha[cell] += coeff_a_vert[vert];
      _coef_beta[cell] += coeff_b_vert[vert];
    }
    _coef_alpha[cell] /= Rectangle::n_vertices;
    _coef_beta[cell] /= Rectangle::n_vertices;
  }
}
