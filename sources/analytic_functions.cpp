#include "analytic_functions.h"
#include "fem/math_functions.h"
#include "fem/auxiliary_functions.h"
#include <iostream>

using namespace fem;

double AnalyticSolution::value(const Point &p, const double t) const
{
  //const double x = p.coord(0);
  //const double y = p.coord(1);
  return 0;
}


RHSFunction::RHSFunction(const Parameters &param)
  : h(param.SOURCE_SUPPORT), // some constant
    f(param.SOURCE_FREQUENCY), // frequency, Hz
    xc(param.SOURCE_CENTER_X), // x-coordinate of the source center
    yc(param.SOURCE_CENTER_Y) // y-coordinate of the source center
{
  if (h < 0)
  {
    const double hx = (param.X_END - param.X_BEG) / param.N_FINE_X; // size of fine grid in x-direction
    const double hy = (param.Y_END - param.Y_BEG) / param.N_FINE_Y; // size of fine grid in y-direction
    require(fabs(hx - hy) < math::FLOAT_NUMBERS_EQUALITY_TOLERANCE, "The cells are not squares, you can't use negative p");
    h = fabs(h) * hx;
  }
}


double RHSFunction::value(const Point &p, const double t) const
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  const double pi = math::PI;
  const double part = pi*pi*f*f*(t - 2./f)*(t - 2./f);
  const double val = 1./(h*h) * exp(-((x - xc)*(x - xc) + (y - yc)*(y - yc)) / (h*h)) * (1. - 2.*part) * exp(-part);
#if defined(DEBUG)
  std::cout << p << ": rhs, part = " << part << " val = " << val << std::endl;
#endif
  return val;
}


double BoundaryFunction::value(const Point &p, const double t) const
{
//  const double x = p.coord(0);
//  const double y = p.coord(1);
  return 0; //analytic_solution(p);
}


double InitialSolution::value(const Point &p, const double t) const
{
  return 0;
}


// ========================================
//
// testing
//
// ========================================

double an_solution_1::value(const Point &p, const double t) const
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  return x + y;
}
double an_rhs_function_1::value(const Point &p, const double t) const
{
//  const double x = p.coord(0);
//  const double y = p.coord(1);
  return 0;
}



double an_solution_2::value(const Point &p, const double t) const
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  return x * y;
}
double an_rhs_function_2::value(const Point &p, const double t) const
{
//  const double x = p.coord(0);
//  const double y = p.coord(1);
  return 0;
}



double an_solution_3::value(const Point &p, const double t) const
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  return x*x + y*y;
}
double an_rhs_function_3::value(const Point &p, const double t) const
{
//  const double x = p.coord(0);
//  const double y = p.coord(1);
  return -4;
}



double an_solution_4::value(const Point &p, const double t) const
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  return sin(x) + sin(y);
}
double an_rhs_function_4::value(const Point &p, const double t) const
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  return sin(x) + sin(y);
}



double an_solution_5::value(const Point &p, const double t) const
{
  const double x = p.coord(0);
//  const double y = p.coord(1);
  return exp(x);
}
double an_rhs_function_5::value(const Point &p, const double t) const
{
  const double x = p.coord(0);
//  const double y = p.coord(1);
  return -exp(x);
}
