#ifndef ANALYTIC_FUNCTIONS_H
#define ANALYTIC_FUNCTIONS_H

#include "fem/point.h"
#include "fem/function.h"
#include "parameters.h"

using namespace fem;

class AnalyticSolution : public Function
{
public:
  double value(const Point &p, const double t = 0) const;
};

class RHSFunction : public Function
{
public:
  RHSFunction(const Parameters &param);
  double value(const Point &p, const double t = 0) const;

private:
  // parameters of Ricker wavelet
  double h;
  double f;
  double xc;
  double yc;
};

class BoundaryFunction : public Function
{
public:
  double value(const Point &p, const double t = 0) const;
};

class InitialSolution : public Function
{
public:
  double value(const Point &p, const double t = 0) const;
};


// ========================================
//
// testing
//
// ========================================

class an_solution_1 : public Function
{
public:
  double value(const Point &p, const double t = 0) const;
};
class an_rhs_function_1 : public Function
{
public:
  double value(const Point &p, const double t = 0) const;
};

class an_solution_2 : public Function
{
public:
  double value(const Point &p, const double t = 0) const;
};
class an_rhs_function_2 : public Function
{
public:
  double value(const Point &p, const double t = 0) const;
};

class an_solution_3 : public Function
{
public:
  double value(const Point &p, const double t = 0) const;
};
class an_rhs_function_3 : public Function
{
public:
  double value(const Point &p, const double t = 0) const;
};

class an_solution_4 : public Function
{
public:
  double value(const Point &p, const double t = 0) const;
};
class an_rhs_function_4 : public Function
{
public:
  double value(const Point &p, const double t = 0) const;
};

class an_solution_5 : public Function
{
public:
  double value(const Point &p, const double t = 0) const;
};
class an_rhs_function_5 : public Function
{
public:
  double value(const Point &p, const double t = 0) const;
};


#endif // ANALYTIC_FUNCTIONS_H
