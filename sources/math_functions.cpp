#include "math_functions.h"

double rel_error(const Vec &vec1, const Vec &vec2)
{
  Vec diff;
  VecDuplicate(vec1, &diff); // allocate the memory for diff
  VecCopy(vec1, diff); // copy vec1 to diff
  VecAXPY(diff, -1, vec2); // diff = vec1 - vec2
  double norm_diff, norm_vec2;
  VecNorm(diff, NORM_2, &norm_diff);
  VecNorm(vec2, NORM_2, &norm_vec2);
  VecDestroy(&diff); // free the memory
  return norm_diff / norm_vec2;
}
