#ifndef ACOUSTIC2D_H
#define ACOUSTIC2D_H

#include "config.h"
#include "fem/fine_mesh.h"
#include "petscvec.h"
#include "petscmat.h"
#include "fem/dof_handler.h"
#include "fem/csr_pattern.h"

class Parameters;


class Acoustic2D
{
public:
            /**
             * Constructor
             * @param param - parameters of the problem
             */
  Acoustic2D(Parameters *param);

            /**
             * Destructor
             */
  ~Acoustic2D();

            /**
             * The main function launching the calculations
             */
  void solve_triangles();
  void solve_rectangles();


private:
            /**
             * Parameters of the problem
             */
  Parameters *_param;

            /**
             * Fine triangular mesh
             */
  FineMesh _fmesh;

            /**
             * Global vector of right hand side
             */
  Vec _global_rhs;

            /**
             * Global mass matrix
             */
  Mat _global_mass_mat;

            /**
             * Global stiffness matrix
             */
  Mat _global_stiff_mat;

            /**
             * Coefficient by time derivative
             */
  std::vector<double> _coef_alpha;

            /**
             * Coefficient under div
             */
  std::vector<double> _coef_beta;

  Acoustic2D(const Acoustic2D&); /** copy constructor */
  Acoustic2D& operator=(const Acoustic2D&); /** copy assignment operator */

            /**
             * Find the list of boundary nodes
             * @param b_nodes - output list of boundary nodes
             */
  //void find_bound_nodes(std::vector<int> &b_nodes) const;

  void solve_explicit_triangles(const DoFHandler &dof_handler, const CSRPattern &csr_pattern);
  void solve_explicit_rectangles(const DoFHandler &dof_handler, const CSRPattern &csr_pattern);
  void solve_crank_nicolson(const DoFHandler &dof_handler, const CSRPattern &csr_pattern);

  void coefficients_initialization();
  void create_bin_layers_file() const;
  void create_ave_layers_file() const;

  void export_coefficients_distribution(const std::string &filename) const;
};


#endif // ACOUSTIC2D_H
