#ifndef BRANCH_BOUND_H
#define BRANCH_BOUND_H

#include "solution.h"
#include "problem.h"
#include "branch_bound/node.h"

#define MAX_N 500
#define MAX_M 500

typedef uint32_t (*solve_fn)(uint32_t n, uint32_t m, uint32_t is_max, const gsl_vector* c, const gsl_matrix* A,
                             const gsl_vector* b, int32_t* B, int32_t* N, solution_t* solution_ptr,
                             uint32_t* iter_n_ptr);

uint32_t solve_relaxation(solve_fn solver, uint32_t is_max, bb_node_t* node_ptr, int32_t* N, solution_t* solution_ptr,
                          uint32_t* iter_n_ptr);

// Choses a non-integer variable to start branching from.
// Returns -2 on error, -1 if the solution contains only
// integers or the index of the first non-integer
// variable on success
int32_t select_branch_var(const var_arr_t* var_arr_ptr, const solution_t* current_sol_ptr);

// Branch and bound method on linear problem p
uint32_t branch_and_bound(problem_t* problem_ptr, solution_t* solution_ptr);

#endif