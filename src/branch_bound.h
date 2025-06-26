#ifndef BRANCH_BOUND_H
#define BRANCH_BOUND_H

#include "problem.h"

// Choses a non-integer variable to start branching from.
// Returns -2 on error, -1 if the solution contains only
// integers or the index of the first non-integer
// variable on success
int32_t problem_select_branch_var(const solution_t s, const uint32_t* integers);

// Creates a duplicate of p with an added constraint
problem_t problem_branch(const problem_t p, int32_t branch_var_index, double bound, char direction);

// Branch and bound method on linear problem p
solution_t branch_and_bound(const problem_t p);

#endif