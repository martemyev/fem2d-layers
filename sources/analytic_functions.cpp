#include "analytic_functions.h"
#include <cmath>

double analytic_solution(const Point &p)
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  return 0; //x + y;
}


double rhs_function(const Point &p, double t, const Parameters &param)
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  const double h = param.SOURCE_SUPPORT; // some constant
  const double f = param.SOURCE_FREQUENCY; // frequency, Hz
  const double xc = param.SOURCE_CENTER_X; // x-coordinate of the source center
  const double yc = param.SOURCE_CENTER_Y; // y-coordinate of the source center
  const double pi = 3.141592654;
  const double part = pi*pi*f*f*(t - 2./f)*(t - 2./f);
  return h*h * exp(-h*h*((x - xc)*(x - xc) + (y - yc)*(y - yc))) * (1. - 2.*part) * exp(-part);
}


double boundary_function(const Point &p, double t)
{
//  const double x = p.coord(0);
//  const double y = p.coord(1);
  return 0; //analytic_solution(p);
}


double init_solution(const Point &p, double t)
{
  return 0;
}


// ========================================
//
// testing
//
// ========================================

double an_solution_1(const Point &p, double t)
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  return x + y;
}
double an_rhs_function_1(const Point &p, double t, const Parameters &param)
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  return 0;
}



double an_solution_2(const Point &p, double t)
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  return x * y;
}
double an_rhs_function_2(const Point &p, double t, const Parameters &param)
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  return 0;
}



double an_solution_3(const Point &p, double t)
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  return x*x + y*y;
}
double an_rhs_function_3(const Point &p, double t, const Parameters &param)
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  return -4;
}



double an_solution_4(const Point &p, double t)
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  return sin(x) + sin(y);
}
double an_rhs_function_4(const Point &p, double t, const Parameters &param)
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  return sin(x) + sin(y);
}



double an_solution_5(const Point &p, double t)
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  return exp(x);
}
double an_rhs_function_5(const Point &p, double t, const Parameters &param)
{
  const double x = p.coord(0);
  const double y = p.coord(1);
  return -exp(x);
}
