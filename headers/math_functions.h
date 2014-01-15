#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

#include "petscvec.h"


/**
 * Relative error in 2-norm between two vectors
 * rel_error = norm2(vec1-vec2) / norm2(vec2)
 */
double rel_error(const Vec &vec1, const Vec &vec2);



#endif // MATH_FUNCTIONS_H
