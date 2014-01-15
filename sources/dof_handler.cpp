#include "dof_handler.h"
#include "fine_mesh.h"
#include "parameters.h"
#include "auxiliary_functions.h"


DoFHandler::DoFHandler(FineMesh *fmesh, const Parameters &param)
  : _fmesh(fmesh)
{
  distribute(param); // distribute dofs
}



DoFHandler::~DoFHandler()
{ }



void DoFHandler::distribute(const Parameters &param)
{
  switch (param.FE_ORDER)
  {
  case 1:
    distribute_first(param);
    return;
  default:
    require(false, "Incorrect order of fe basis functions (" + d2s(param.FE_ORDER) + ")");
  }
}



void DoFHandler::distribute_first(const Parameters &param)
{
  require(param.FE_ORDER == 1, "This is not implemented for fe order other than 1");

  _dofs.resize(_fmesh->n_vertices());
  for (int ver = 0; ver < _fmesh->n_vertices(); ++ver)
    _dofs[ver] = _fmesh->vertex(ver);

  for (int tr = 0; tr < _fmesh->n_triangles(); ++tr)
  {
    Triangle *triangle = _fmesh->triangle_orig(tr);
    triangle->n_dofs(Triangle::n_dofs_first); // allocate the memory for degrees of freedom
    for (int ver = 0; ver < Triangle::n_vertices; ++ver)
    {
      // set the numbers of degrees of freedom.
      // the number of the degree of freedom is the same as the number of the corresponding vertex
      triangle->dof(ver, triangle->vertex(ver));
    }
  }
}



unsigned int DoFHandler::n_dofs() const
{
  return _dofs.size();
}



Point DoFHandler::dof(unsigned int number) const
{
  expect(number >= 0 && number < _dofs.size(), "Incorrect input parameter");
  return _dofs[number];
}



const FineMesh* DoFHandler::fmesh() const
{
  return _fmesh;
}



