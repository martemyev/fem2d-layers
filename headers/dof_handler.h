#ifndef DOF_HANDLER_H
#define DOF_HANDLER_H

#include "point.h"
#include <vector>


class FineMesh;
class Parameters;

class DoFHandler
{
public:
            /**
             * Constructor
             */
  DoFHandler(FineMesh *fmesh, const Parameters &param);

            /**
             * Destructor
             */
  ~DoFHandler();

            /**
             * Get the number of the degrees of freedom
             */
  unsigned int n_dofs() const;

            /**
             * Get the degree of freedom (its copy)
             * @param number - the serial number of the dof
             */
  Point dof(unsigned int number) const;

            /**
             * Const pointer to the fine mesh for access to mesh's functions
             */
  const FineMesh* fmesh() const;


private:
            /**
             * Finite triangular mesh which the dof handler is connected with
             */
  FineMesh *_fmesh;

            /**
             * Degrees of freedom
             */
  std::vector<Point> _dofs;

            /**
             * Copy constructor is private to prevent copying of dof handlers
             */
  DoFHandler(const DoFHandler&);

            /**
             * Copy assignment operator is also private to prevent copying
             */
  DoFHandler& operator=(const DoFHandler&);

            /**
             * Distribute degrees of freedom according to fine mesh elements (connectivity of vertices),
             * order of basis functions (from param object), and type of finite element (CG, DG, etc).
             * We consider only CG case right now though.
             */
  void distribute(const Parameters &param);

            /**
             * Distribute degrees of freedom in case of first order basis functions
             */
  void distribute_first(const Parameters &param);
};


#endif // DOF_HANDLER_H
