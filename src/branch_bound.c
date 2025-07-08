// #include "branch_bound.h"
// #include "simplex.h"
// #include "pstack.h"
// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
// #include <gsl/gsl_vector.h>
// #include <gsl/gsl_matrix.h>

// // Choses a non-integer variable to start branching from.
// // Returns -2 on error, -1 if the solution contains only
// // integers or the index of the first non-integer
// // variable on success
// int32_t problem_select_branch_var(const varr_t* varr, uint32_t variables_num, const solution_t* current_solution) {
//     if (!varr || !current_solution) {
//         fprintf(stderr, "Some arguments are NULL in problem_select_branch_var\n");
//         return -2;
//     }

//     for (uint32_t i = 0; i < variables_num; i++) {
//         if (variable_is_integer(varr_get(varr, i)) && !solution_var_is_integer(current_solution, i)) {
//             return (int32_t)i;
//         }
//     }

//     return -1;
// }

// // Creates a duplicate of p with an added bound on the variable
// // with index branch_var_index.
// // Direction can either be 'L' (Lower bound) or 'U' (Upper bound)
// problem_t* problem_branch(const problem_t* p, int32_t branch_var_index, double bound, char direction) {
//     if (!p || (direction != 'U' && direction != 'L')) {
//         fprintf(stderr, "NULL problem or unknown direction in problem_branch\n");
//         return NULL;
//     }

//     uint32_t n = problem_n(p);
//     uint32_t m = problem_m(p);
//     uint32_t is_max = problem_is_max(p);
//     const gsl_vector* c = problem_c(p);
//     const gsl_vector* b = problem_b(p);
//     const gsl_matrix* A = problem_A(p);
//     const int32_t* basis = problem_basis(p);
//     const varr_t* varr = problem_varr(p);

//     if (!c || !b || !A || !basis || !varr) {
//         fprintf(stderr, "Error in getting resources from problem p in problem_branch\n");
//         return NULL;
//     }

//     uint32_t n2 = n + 1;  // +1 for the new constraint
//     uint32_t m2 = m + 1;  // +1 for the new slack/surplus variable

//     gsl_vector* c2 = NULL;
//     gsl_matrix* A2 = NULL;
//     gsl_vector* b2 = NULL;
//     variable_t** variables2 = NULL;
//     int32_t* basis2 = NULL;

//     c2 = gsl_vector_alloc(m2);
//     if (!c2) {
//         goto fail;
//     }

//     A2 = gsl_matrix_alloc(n2, m2);
//     if (!A2) {
//         goto fail;
//     }

//     b2 = gsl_vector_alloc(n2);
//     if (!b2) {
//         goto fail;
//     }

//     uint32_t vars_set = 0;
//     variables2 = (variable_t**)malloc(sizeof(variable_t*) * m2);
//     if (!variables2) {
//         goto fail;
//     }

//     // Direction == 'U' => x[branch_var_index] <= floor(bound)
//     // Direction == 'L' => x[branch_var_index] >= ceil(bound)
//     bound = direction == 'U' ? floor(bound) : ceil(bound);

//     // Copy c and variables2 adding new values
//     for (uint32_t j = 0; j < m2; j++) {
//         if (j < m) {
//             // Copy old c value
//             gsl_vector_set(c2, j, gsl_vector_get(c, j));

//             // Copy old variable
//             variables2[j] = variable_duplicate(varr_get(varr, j));
//         } else {
//             // Set 0 for new slack/surplus variable
//             gsl_vector_set(c2, j, 0.0);

//             // Add new slack/surplus variable
//             variables2[j] = variable_new_real_positive(10e9);
//         }

//         vars_set++;
//         if (!variables2[j]) {
//             goto fail;
//         }
//     }

//     // Copy b and A adding new values
//     for (uint32_t i = 0; i < n2; i++) {
//         if (i < n) {
//             // Copy old b value
//             gsl_vector_set(b2, i, gsl_vector_get(b, i));
//         } else {
//             // Add bound value
//             gsl_vector_set(b2, n, bound);
//         }

//         for (uint32_t j = 0; j < m2; j++) {
//             if (i < n && j == m) {
//                 // Set 0 for the new slack/surplus variable
//                 // in every old constraint
//                 gsl_matrix_set(A2, i, j, 0.0);
//             } else if (i == n) {
//                 // Set 0 for every variable in the new
//                 // constraint except for the branch
//                 // variable and the new slack/surplus variable
//                 if (j == (uint32_t)branch_var_index) {
//                     gsl_matrix_set(A2, i, j, 1.0);
//                 } else if (j < m) {
//                     gsl_matrix_set(A2, i, j, 0.0);
//                 } else if (j == m) {
//                     // direction == 'U' => slack variable (+)
//                     // direction == 'L' => surplus variable (-)
//                     gsl_matrix_set(A2, i, j, direction == 'U' ? 1.0 : -1.0);
//                 }
//             } else {
//                 // Copy old A value
//                 gsl_matrix_set(A2, i, j, gsl_matrix_get(A, i, j));
//             }
//         }
//     }

//     // Check if original base is still feasible
//     uint32_t og_base_feasible = 1;
//     for (uint32_t i = 0; og_base_feasible && i < n; i++) {
//         double bi = gsl_vector_get(b2, i);
//         if (bi < -1e-8) {
//             og_base_feasible = 0;
//         }
//     }

//     if (og_base_feasible && basis) {
//         // Copy the base adding the new slack/surplus variable
//         basis2 = (int32_t*)malloc(sizeof(int32_t) * n2);
//         if (!basis2) {
//             goto fail;
//         }

//         for (uint32_t i = 0; i < n; i++) {
//             basis2[i] = basis[i];
//         }

//         // New basis index
//         basis2[n] = m2 - 1;
//     } else {
//         // Find a new feasible base
//         uint32_t pI_iters = 0;
//         basis2 = simplex_phaseI(n2, m2, A2, b2, variables2, &pI_iters);
//         if (!basis2) {
//             fprintf(stderr, "Failed to find basis with Phase1 in problem_branch\n");
//             goto fail;
//         }
//     }

//     return problem_new(n2, m2, is_max, c2, A2, b2, basis2, problem_pI_iter(p), variables2);

// fail:
//     gsl_vector_free(c2);
//     gsl_matrix_free(A2);
//     gsl_vector_free(b2);
//     if (variables2) {
//         for (uint32_t i = 0; i < vars_set; i++) {
//             variable_free(&variables2[i]);
//         }
//         free(variables2);
//     }
//     free(basis2);
//     return NULL;
// }

// // Branch and bound method on linear problem p
// solution_t* branch_and_bound(const problem_t* p) {
//     if (!p) {
//         return NULL;
//     }

//     pstack_t* stack = pstack_new();
//     if (!pstack_push(stack, problem_duplicate(p))) {
//         pstack_free(&stack);
//         return NULL;
//     }

//     solution_t* best = NULL;
//     const varr_t* varr = problem_varr(p);
//     uint32_t m = problem_m(p);

//     while (!pstack_empty(stack)) {
//         problem_t* current_problem = pstack_pop(stack);
//         if (!current_problem) {
//             continue;
//         }

//         solution_t* current_solution = simplex_phaseII(current_problem);

//         if (!current_solution || solution_is_unbounded(current_solution)) {
//             solution_free(&current_solution);
//             problem_free(&current_problem);
//             continue;
//         }

//         int32_t branch_var = problem_select_branch_var(varr, m, current_solution);
//         if (branch_var == -2) {
//             solution_free(&current_solution);
//             problem_free(&current_problem);
//             continue;
//         }

//         // Only integer variables
//         if (branch_var == -1) {
//             if (!best) {
//                 best = current_solution;
//             } else if (problem_is_max(current_problem) ? solution_z(current_solution) > solution_z(best)
//                                                        : solution_z(current_solution) < solution_z(best)) {
//                 solution_free(&best);
//                 best = current_solution;
//             } else {
//                 solution_free(&current_solution);
//             }

//             problem_free(&current_problem);
//             continue;
//         }

//         // Branch on non-integer variable
//         double bound = gsl_vector_get(solution_x(current_solution), branch_var);
//         problem_t* left = problem_branch(current_problem, branch_var, bound, 'U');
//         problem_t* right = problem_branch(current_problem, branch_var, bound, 'L');

//         if (left && !pstack_push(stack, left)) {
//             problem_free(&left);
//         }

//         if (right && !pstack_push(stack, right)) {
//             problem_free(&right);
//         }

//         solution_free(&current_solution);
//         problem_free(&current_problem);
//     }

//     pstack_free(&stack);

//     return best;
// }