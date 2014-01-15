#ifndef RESULT_H
#define RESULT_H

#include "petscvec.h"
#include <string>

class DoFHandler;


class Result
{
public:
            /**
             * Constructor
             * @param dof_handler - pointer to the handler of degrees of freedom
             */
  Result(const DoFHandler *dof_handler);

            /**
             * Destructor
             */
  ~Result();

            /**
             * Write the results to file in vtk(vtu) format to work then in Paraview
             * @param filename - the name of output file
             */
  void write_vtu(const std::string &filename, const Vec &solution, const Vec &exact_solution = 0) const;


private:
            /**
             * Constant pointer to dof_handler
             */
  const DoFHandler *_dof_handler;

            /**
             * Copy constructor
             */
  Result(const Result&);

            /**
             * Copy assignment operator
             */
  Result& operator=(const Result&);
};


#endif // RESULT_H
