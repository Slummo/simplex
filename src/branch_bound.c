#include "branch_bound.h"
#include "simplex.h"
#include "pstack.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

// Choses a non-integer variable to start branching from.
// Returns -2 on error, -1 if the solution contains only
// integers or the index of the first non-integer
// variable on success
int32_t problem_select_branch_var(const solution_t s) {
    if (!s) {
        fprintf(stderr, "s is NULL in problem_select_branch_var\n");
        return -2;
    }

    for (uint32_t i = 0; i < solution_x(s)->size; i++) {
        if (!solution_var_is_integer(s, i)) {
            return (int32_t)i;
        }
    }

    return -1;
}

// Creates a duplicate of p with an added bound on the variable
// with index branch_var_index.
// Direction can either be 'L' (Lower bound) or 'U' (Upper bound)
problem_t problem_branch(const problem_t p, int32_t branch_var_index, double bound, char direction) {
    if (!p || (direction != 'U' && direction != 'L')) {
        fprintf(stderr, "Error in problem_branch\n");
        return NULL;
    }

    uint32_t n = problem_n(p);
    uint32_t m = problem_m(p);
    const gsl_vector* c = problem_c(p);
    const gsl_vector* b = problem_b(p);
    const gsl_matrix* A = problem_A(p);
    const int32_t* basis = problem_basis(p);

    uint32_t n2 = n + 1;  // +1 for the new constraint
    uint32_t m2 = m + 1;  // +1 for the new slack/surplus variable

    gsl_vector* c2 = gsl_vector_alloc(m2);
    if (!c2) {
        return NULL;
    }

    gsl_vector* b2 = gsl_vector_alloc(n2);
    if (!b2) {
        gsl_vector_free(c2);
        return NULL;
    }

    gsl_matrix* A2 = gsl_matrix_alloc(n2, m2);
    if (!A2) {
        gsl_vector_free(c2);
        gsl_vector_free(b2);
        return NULL;
    }

    // Direction == 'U' => x[branch_var_index] <= floor(bound)
    // Direction == 'L' => x[branch_var_index] >= ceil(bound)
    bound = direction == 'U' ? floor(bound) : ceil(bound);

    // Copy c adding new value
    for (uint32_t j = 0; j < m2; j++) {
        if (j < m) {
            // Copy old c value
            gsl_vector_set(c2, j, gsl_vector_get(c, j));
        } else {
            // Set 0 for new slack/surplus variable
            gsl_vector_set(c2, j, 0.0);
        }
    }

    // Copy b and A adding new values
    for (uint32_t i = 0; i < n2; i++) {
        if (i < n) {
            // Copy old b value
            gsl_vector_set(b2, i, gsl_vector_get(b, i));
        } else {
            // Add bound value
            gsl_vector_set(b2, n, bound);
        }

        for (uint32_t j = 0; j < m2; j++) {
            if (i < n && j == m) {
                // Set 0 for the new slack/surplus variable
                // in every old constraint
                gsl_matrix_set(A2, i, j, 0.0);
            } else if (i == n) {
                // Set 0 for every variable in the new
                // constraint except for the branch
                // variable and the new slack/surplus variable
                if (j == (uint32_t)branch_var_index) {
                    gsl_matrix_set(A2, i, j, 1.0);
                } else if (j < m) {
                    gsl_matrix_set(A2, i, j, 0.0);
                } else if (j == m) {
                    // direction == 'U' => slack variable (+)
                    // direction == 'L' => surplus variable (-)
                    gsl_matrix_set(A2, i, j, direction == 'U' ? 1.0 : -1.0);
                }
            } else {
                // Copy old A value
                gsl_matrix_set(A2, i, j, gsl_matrix_get(A, i, j));
            }
        }
    }

    int32_t* basis2 = NULL;

    // Check if original base is still feasible
    uint32_t og_base_feasible = 1;
    for (uint32_t i = 0; i < n; i++) {
        double bi = gsl_vector_get(b2, i);
        if (bi < -1e-8) {
            og_base_feasible = 0;
            break;
        }
    }

    if (og_base_feasible && basis) {
        // Copy the base adding the new slack/surplus variable
        basis2 = (int32_t*)malloc(sizeof(int32_t) * n2);
        if (!basis2) {
            gsl_vector_free(c2);
            gsl_vector_free(b2);
            gsl_matrix_free(A2);
            return NULL;
        }

        for (uint32_t i = 0; i < n; i++) {
            basis2[i] = basis[i];
        }

        // New basis index
        basis2[n] = m2 - 1;
    } else {
        // Find a new feasible base
        uint32_t pI_iters = 0;
        basis2 = simplex_phaseI(n2, m2, A2, b2, &pI_iters);
        if (!basis2) {
            fprintf(stderr, "Failed to find basis with Phase1 in problem_branch\n");
            gsl_vector_free(c2);
            gsl_vector_free(b2);
            gsl_matrix_free(A2);
            return NULL;
        }
    }

    return problem_new(n2, m2, problem_is_max(p), c2, A2, b2, basis2, problem_pI_iter(p));
}

// Branch and bound method on linear problem p
solution_t branch_and_bound(const problem_t p) {
    if (!p) {
        return NULL;
    }

    pstack_t stack = pstack_new();
    if (!stack) {
        return NULL;
    }

    if (!pstack_push(stack, p)) {
        pstack_free(&stack);
        return NULL;
    }

    solution_t best = NULL;

    while (!pstack_empty(stack)) {
        problem_t current = pstack_pop(stack);
        if (!current) {
            continue;
        }

        solution_t s = simplex_phaseII(current);

        if (!s || solution_is_unbounded(s)) {
            solution_free(&s);
            problem_free(&current);
            continue;
        }

        int32_t branch_var = problem_select_branch_var(s);
        if (branch_var == -2) {
            solution_free(&s);
            problem_free(&current);
            continue;
        }

        // Only integer variables
        if (branch_var == -1) {
            if (!best) {
                best = s;
            } else if (problem_is_max(current) ? solution_z(s) > solution_z(best) : solution_z(s) < solution_z(best)) {
                solution_free(&best);
                best = s;
            } else {
                solution_free(&s);
            }

            problem_free(&current);
            continue;
        }

        // Branch on non-integer variable
        double bound = gsl_vector_get(solution_x(s), branch_var);
        problem_t left = problem_branch(current, branch_var, bound, 'U');
        problem_t right = problem_branch(current, branch_var, bound, 'L');

        if (left) {
            pstack_push(stack, left);
        }

        if (right) {
            pstack_push(stack, right);
        }

        solution_free(&s);
        problem_free(&current);
    }

    pstack_free(&stack);

    return best;
}