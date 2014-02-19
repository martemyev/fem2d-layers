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
  void solve();


private:
            /**
             * parameters of the problem
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

  Acoustic2D(const Acoustic2D&); /** copy constructor */
  Acoustic2D& operator=(const Acoustic2D&); /** copy assignment operator */

            /**
             * Find the list of boundary nodes
             * @param b_nodes - output list of boundary nodes
             */
  void find_bound_nodes(std::vector<int> &b_nodes) const;

  void solve_explicit(const DoFHandler &dof_handler, const CSRPattern &csr_pattern);
  void solve_crank_nicolson(const DoFHandler &dof_handler, const CSRPattern &csr_pattern);
};


#endif // ACOUSTIC2D_H
