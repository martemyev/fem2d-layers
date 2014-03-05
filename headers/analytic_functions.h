#ifndef ANALYTIC_FUNCTIONS_H
#define ANALYTIC_FUNCTIONS_H

#include "fem/point.h"
#include "parameters.h"

double analytic_solution(const Point &p);
double rhs_function(const Point &p, double time);
double boundary_function(const Point &p, double t);
double init_solution(const Point &p, double t);


// ========================================
//
// testing
//
// ========================================

double an_solution_1(const Point &p, double t = 0.);
double an_rhs_function_1(const Point &p, double t);

double an_solution_2(const Point &p, double t = 0.);
double an_rhs_function_2(const Point &p, double t);

double an_solution_3(const Point &p, double t = 0.);
double an_rhs_function_3(const Point &p, double t);

double an_solution_4(const Point &p, double t = 0.);
double an_rhs_function_4(const Point &p, double t);

double an_solution_5(const Point &p, double t = 0.);
double an_rhs_function_5(const Point &p, double t);

#endif // ANALYTIC_FUNCTIONS_H
