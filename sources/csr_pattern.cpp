#include "csr_pattern.h"
#include "dof_handler.h"
#include "triangle.h"
#include "fine_mesh.h"
#include "auxiliary_functions.h"
#include <set>


CSRPattern::CSRPattern(const DoFHandler &dof_handler)
{
  make_sparse_format(dof_handler);
}



CSRPattern::~CSRPattern()
{
  _row.clear();
  _col.clear();
}



void CSRPattern::make_sparse_format(const DoFHandler &dof_handler)
{
  // the number of rows of the matrix (and its order) of connectivity between degrees of freedom
  _order = dof_handler.n_dofs();

  std::set<unsigned int> *connect = new std::set<unsigned int>[_order];

  // pass through all triangles and all dofs on them
  for (int cell = 0; cell < dof_handler.fmesh()->n_triangles(); ++cell)
  {
    Triangle triangle = dof_handler.fmesh()->triangle(cell);

    for (int di = 0; di < triangle.n_dofs(); ++di)
    {
      const unsigned int dof_i = triangle.dof(di); // the number of the first degree of freedom
      for (int dj = 0; dj < triangle.n_dofs(); ++dj)
      {
        const unsigned int dof_j = triangle.dof(dj); // the number of the second degree of freedom
        // insert the values in the corresponding places
        connect[dof_i].insert(dof_j);
        connect[dof_j].insert(dof_i);
      }
    }
  }

  // initialization of the pattern
  _row.resize(_order + 1);
  _row[0] = 0;
  for (int i = 0; i < _order; ++i)
    _row[i + 1] = _row[i] + connect[i].size();

  _col.resize(_row[_order]);
  int k = 0;
  for (int i = 0; i < _order; ++i)
  {
    for (std::set<unsigned int>::const_iterator iter = connect[i].begin();
         iter != connect[i].end();
         ++iter)
    {
      _col[k] = *iter;
      ++k;
    }
  }

  // free the memory
  for (int i = 0; i < _order; ++i)
    connect[i].clear();
  delete[] connect;
}



unsigned int CSRPattern::order() const
{
  return _order;
}



unsigned int CSRPattern::row(unsigned int number) const
{
  expect(number >= 0 && number < _row.size(), "Incorrect input data");
  return _row[number];
}



unsigned int CSRPattern::col(unsigned int number) const
{
  expect(number >= 0 && number < _col.size(), "Incorrect input data");
  return _col[number];
}



const int* CSRPattern::nnz() const
{
  int *nnz = new int[_order];
  for (int i = 0; i < _order; ++i)
    nnz[i] = _row[i + 1] - _row[i];
  return nnz;
}
