#include "analytic_functions.h"
#include "fem/math_functions.h"

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
{ }
double RHSFunction::value(const Point &p, const double t) const
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  const double pi = 3.141592654;
  const double part = pi*pi*f*f*(t - 2./f)*(t - 2./f);
  return h*h * exp(-h*h*((x - xc)*(x - xc) + (y - yc)*(y - yc))) * (1. - 2.*part) * exp(-part);
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
