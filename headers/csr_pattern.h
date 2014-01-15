#ifndef SPARSE_FORMAT_H
#define SPARSE_FORMAT_H

#include <vector>

class DoFHandler;

/**
 * Compressed sparse row (CSR) pattern
 */
class CSRPattern
{
public:
            /**
             * Constructor
             */
  CSRPattern(const DoFHandler &dof_handler);

            /**
             * Destructor
             */
  ~CSRPattern();

            /**
             * Get matrix order
             */
  unsigned int order() const;

            /**
             * Get an element from row array
             * @param number - the serial number of the element
             */
  unsigned int row(unsigned int number) const;

            /**
             * Get an element from col array
             * @param number - the serial number of the element
             */
  unsigned int col(unsigned int number) const;

  const int* nnz() const;


private:
            /**
             * Matrix order (matrix is square).
             */
  unsigned int _order;
            /**
             * The number of nonzero elements in each row. Starts from 0.
             * Its dimension = the number of matrix rows + 1.
             * Its last element means the number of nonzero elements in the matrix.
             */
  std::vector<unsigned int> _row;

            /**
             * The numbers of columns where nonzero elements take place.
             * Its dimension is equal to the amount of nonzero elements in the matrix.
             */
  std::vector<unsigned int> _col;

            /**
             * Copy constructor. It's private, since we don't want to copy CSRPattern objects
             */
  CSRPattern(const CSRPattern&);

            /**
             * The same is valid as for copy constructor
             */
  CSRPattern& operator=(const CSRPattern&);

            /**
             * Make sparse format
             * @param dof_handler - handler of degrees of freedom
             */
  void make_sparse_format(const DoFHandler &dof_handler);
};

#endif // SPARSE_FORMAT_H
