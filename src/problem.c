#include "problem.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_linalg.h>

/// Return the smaller of a and b
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

struct problem {
    uint32_t n;          // Number of constraints
    uint32_t m;          // Number of variables
    uint32_t is_max;     // Boolean value to know if its a maximization problem
    gsl_vector* c;       // Reduced costs (m)
    gsl_matrix* A;       // Constraints matrix (n x m)
    gsl_vector* b;       // RHS (n)
    uint32_t* basis;     // Indices of basic variables (size n)
    uint32_t* nonbasis;  // Indices of nonbasic variables (size m-n)
};

// Find problem basis indices with Phase 1 method
uint32_t* problem_find_basis(uint32_t n, uint32_t m, const gsl_matrix* A, const gsl_vector* b, uint32_t* iter_ptr) {
    // One artificial variable per constraint
    gsl_vector* c = gsl_vector_alloc(m + n);
    if (!c) {
        fprintf(stderr, "Failed to alloc c for phase 1\n");
        return NULL;
    }

    // Set vector c
    for (uint32_t i = 0; i < m + n; i++) {
        gsl_vector_set(c, i, i < m ? 0.0 : -1.0);
    }

    uint32_t* artificial_indices = (uint32_t*)malloc(sizeof(uint32_t) * n);
    if (!artificial_indices) {
        fprintf(stderr, "Failed to allocate array for artificial indices\n");
        gsl_vector_free(c);
        return NULL;
    }

    // Set array
    for (uint32_t i = 0; i < n; i++) {
        artificial_indices[i] = m + i;
    }

    // Augmented matrix with m + n columns (one column for each artificial variable)
    gsl_matrix* A2 = gsl_matrix_alloc(n, m + n);
    if (!A2) {
        fprintf(stderr, "Failed to alloc A2\n");
        gsl_vector_free(c);
        free(artificial_indices);
        return NULL;
    }

    // Copy first m values for each row,
    // then add values for artificial variables
    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t j = 0; j < m + n; j++) {
            if (j < m) {
                gsl_matrix_set(A2, i, j, gsl_matrix_get(A, i, j));
            } else if (j == m + i) {
                // Artificial variable in column m + i
                gsl_matrix_set(A2, i, j, 1.0);
            } else {
                gsl_matrix_set(A2, i, j, 0.0);
            }
        }
    }

    gsl_vector* b2 = vector_duplicate(b);
    if (!b2) {
        gsl_vector_free(c);
        free(artificial_indices);
        gsl_matrix_free(A2);
        return NULL;
    }

    problem_t phase1 = problem_new(n, m + n, 1, c, A2, b2, artificial_indices);
    if (!phase1) {
        fprintf(stderr, "Failed to create phase1 problem\n");
        return NULL;
    }

    problem_print(phase1, "Phase1");

    solution_t s = solve(phase1);
    if (!s) {
        fprintf(stderr, "Failed to create solution for phase1 problem\n");
        problem_free(&phase1);
        return NULL;
    }

    uint32_t* basis = (uint32_t*)malloc(sizeof(uint32_t) * n);
    if (!basis) {
        fprintf(stderr, "Failed to allocate array for basis indices\n");
        problem_free(&phase1);
        solution_free(&s);
        return NULL;
    }

    // No feasible base for original problem
    if (solution_optimal_value(s) < 0) {
        free(basis);
        basis = NULL;
    } else {
        gsl_vector* x = solution_optimal_vec(s);
        uint32_t* final_basis = solution_basis(s);
        for (uint32_t i = 0; i < n; i++) {
            // If the variable in basis is artificial (index >= m), make sure its value is 0
            if (final_basis[i] >= m && gsl_vector_get(x, final_basis[i]) > 0) {
                // Infeasible
                free(basis);
                basis = NULL;
                break;
            }

            basis[i] = final_basis[i];
        }
    }

    *iter_ptr = solution_iterations(s);

    problem_free(&phase1);
    solution_free(&s);

    return basis;
}

// Frees arguments on error
problem_t problem_new(uint32_t n, uint32_t m, uint32_t is_max, gsl_vector* c, gsl_matrix* A, gsl_vector* b,
                      uint32_t* basis) {
    problem_t p = (problem_t)malloc(sizeof(_problem));
    if (!p) {
        if (c)
            gsl_vector_free(c);
        if (A)
            gsl_matrix_free(A);
        if (b)
            gsl_vector_free(b);
        if (basis)
            free(basis);
        return NULL;
    }

    if (!is_max) {
        // Multiply by -1 to have the objective function
        // in standard form
        gsl_vector_scale(c, -1.0);
    }

    p->n = n;
    p->m = m;
    p->is_max = is_max;
    p->c = c;
    p->A = A;
    p->b = b;
    p->basis = basis;

    p->nonbasis = (uint32_t*)malloc(sizeof(uint32_t) * (m - n));
    if (!p->nonbasis) {
        problem_free(&p);
        return NULL;
    }

    // Fill non-basis from remaining indices
    uint32_t* used = calloc(m, sizeof(uint32_t));
    for (uint32_t i = 0; i < n; i++) {
        if (used[basis[i]]) {
            fprintf(stderr, "Duplicate basis index %u\n", basis[i]);
            problem_free(&p);
            return NULL;
        }
        used[basis[i]] = 1;
    }

    // Fill non-basis with those not used
    uint32_t ni = 0;
    for (uint32_t i = 0; i < m; i++) {
        if (!used[i]) {
            p->nonbasis[ni++] = i;
        }
    }
    free(used);

    return p;
}

#define TERM_WIDTH 8

void print_coefficient(double coef, uint32_t var_idx, int is_first) {
    char buf[TERM_WIDTH + 1];

    if (fabs(coef) < 1e-9) {
        // Print spaces for alignment
        snprintf(buf, sizeof(buf), "%*s", TERM_WIDTH, "");
    } else {
        char sign = coef < 0 ? '-' : (is_first ? ' ' : '+');
        double abs_val = fabs(coef);

        // Check if it's 1 or -1
        if (fabs(abs_val - 1.0) < 1e-9) {
            snprintf(buf, sizeof(buf), "%c x%-2u", sign, var_idx + 1);  // omit the 1
        } else if (fabs(abs_val - (int32_t)abs_val) < 1e-9) {
            snprintf(buf, sizeof(buf), "%c%2dx%-2u", sign, (int32_t)abs_val, var_idx + 1);
        } else {
            snprintf(buf, sizeof(buf), "%c%.1lfx%-2u", sign, abs_val, var_idx + 1);
        }
    }

    printf("%*s", TERM_WIDTH, buf);
}

void problem_print(const problem_t p, const char* name) {
    if (!p)
        return;

    printf("================== %s ==================\n", name);
    printf("%s z = ", p->is_max ? "max" : "min");

    // Objective function
    for (uint32_t i = 0; i < p->m; i++) {
        double ci = gsl_vector_get(p->c, i);
        print_coefficient(ci, i, i == 0);
    }
    printf("\n\nconstraints:\n");

    // Constraints
    for (uint32_t i = 0; i < p->n; i++) {
        printf("\t");
        for (uint32_t j = 0; j < p->m; j++) {
            double aij = gsl_matrix_get(p->A, i, j);
            print_coefficient(aij, j, j == 0);
        }

        double bi = gsl_vector_get(p->b, i);
        printf(" = ");
        if (fabs(bi - (int32_t)bi) < 1e-9) {
            printf("%d", (int32_t)bi);
        } else {
            printf("%.3lf", bi);
        }
        printf("\n");
    }

    printf("\n");
}

void problem_free(problem_t* pp) {
    if (!pp || !*pp) {
        return;
    }

    if ((*pp)->c)
        gsl_vector_free((*pp)->c);
    if ((*pp)->A)
        gsl_matrix_free((*pp)->A);
    if ((*pp)->b)
        gsl_vector_free((*pp)->b);
    if ((*pp)->basis)
        free((*pp)->basis);
    if ((*pp)->nonbasis)
        free((*pp)->nonbasis);
    free(*pp);
    *pp = NULL;
}

solution_t solve(problem_t problem) {
    if (!problem) {
        return NULL;
    }

    uint32_t n = problem->n;
    uint32_t m = problem->m;
    gsl_matrix* A = problem->A;
    gsl_vector* b = problem->b;
    gsl_vector* c = problem->c;
    uint32_t* B = problem->basis;
    uint32_t* N = problem->nonbasis;

    gsl_matrix* Ab = gsl_matrix_alloc(n, n);
    gsl_matrix* Ab_inv = NULL;
    gsl_vector* xB = gsl_vector_alloc(n);
    gsl_vector* cb = gsl_vector_alloc(n);
    gsl_vector* cn = gsl_vector_alloc(m - n);
    gsl_vector* r = gsl_vector_alloc(m - n);

    if (!Ab || !xB || !cb || !cn || !r) {
        if (Ab)
            gsl_matrix_free(Ab);
        if (xB)
            gsl_vector_free(xB);
        if (cb)
            gsl_vector_free(cb);
        if (cn)
            gsl_vector_free(cn);
        if (r)
            gsl_vector_free(r);
        return NULL;
    }

    uint32_t iter = 0;
    uint32_t unbounded = 0;
    while (1) {
        // Extract Ab matrix
        for (uint32_t i = 0; i < n; i++) {
            for (uint32_t j = 0; j < n; j++) {
                // B indices must be valid column indices
                gsl_matrix_set(Ab, i, j, gsl_matrix_get(A, i, B[j]));
            }
            gsl_vector_set(cb, i, gsl_vector_get(c, B[i]));
        }

        // Compute the inverse of Ab
        Ab_inv = inverse(Ab, n);

        // Compute xB = Ab_inv * b
        gsl_blas_dgemv(CblasNoTrans, 1.0, Ab_inv, b, 0.0, xB);

        // Compute reduced costs for all nonbasis variables
        for (uint32_t i = 0; i < (m - n); i++) {
            uint32_t j = N[i];

            // cn[i] = cj
            gsl_vector_set(cn, i, gsl_vector_get(c, j));

            // aj
            gsl_vector* aj = gsl_vector_alloc(n);

            for (uint32_t k = 0; k < n; k++) {
                gsl_vector_set(aj, k, gsl_matrix_get(A, k, j));
            }

            // Ab_inv * aj
            gsl_vector* Ab_inv_aj = gsl_vector_alloc(n);

            gsl_blas_dgemv(CblasNoTrans, 1.0, Ab_inv, aj, 0.0, Ab_inv_aj);

            // r[i] = rj = cj - cb * Ab_inv * aj
            double dot;
            gsl_blas_ddot(cb, Ab_inv_aj, &dot);
            gsl_vector_set(r, i, gsl_vector_get(cn, i) - dot);

            gsl_vector_free(aj);
            gsl_vector_free(Ab_inv_aj);
        }

        // Choose the entering variable
        int32_t entering = -1;
        for (uint32_t i = 0; i < (m - n); i++) {
            if (gsl_vector_get(r, i) > 0) {
                entering = i;
                break;  // Bland's rule: choose first positive reduced cost
            }
        }

        if (entering == -1) {
            break;  // Optimal
        }

        // Extract column j of A
        uint32_t j = N[entering];
        gsl_vector* aj = gsl_vector_alloc(n);
        for (uint32_t k = 0; k < n; k++) {
            gsl_vector_set(aj, k, gsl_matrix_get(A, k, j));
        }

        // Compute direction vector d = -Ab_inv * aj
        gsl_vector* d = gsl_vector_alloc(n);

        gsl_matrix_scale(Ab_inv, -1.0);
        gsl_blas_dgemv(CblasNoTrans, 1.0, Ab_inv, aj, 0.0, d);

        // Choose leaving variable
        double min_ratio = 1e20;
        int32_t leaving = -1;
        for (uint32_t i = 0; i < n; i++) {
            double di = gsl_vector_get(d, i);
            // Only include variables with negative direction coefficient
            if (di < 0) {
                double ratio = -gsl_vector_get(xB, i) / di;
                if (ratio < min_ratio) {
                    min_ratio = ratio;
                    leaving = i;
                }
            }
        }

        if (leaving == -1) {
            unbounded = 1;
            gsl_vector_free(aj);
            gsl_vector_free(d);
            break;
        }

        // Switch entering and leaving variables
        uint32_t tmp = B[leaving];
        B[leaving] = N[entering];
        N[entering] = tmp;

        gsl_vector_free(aj);
        gsl_vector_free(d);
        iter++;
    }

    solution_t s = solution_new(m, gsl_vector_calloc(m), B, unbounded, iter);

    // Extract optimal solution and value
    if (s && !is_solution_unbounded(s)) {
        gsl_vector* v = solution_optimal_vec(s);
        double z = 0.0;
        for (uint32_t i = 0; i < n; i++) {
            gsl_vector_set(v, B[i], gsl_vector_get(xB, i));
        }
        gsl_blas_ddot(c, v, &z);
        solution_set_optimal_value(s, z);
    }

    gsl_vector_free(cb);
    gsl_vector_free(cn);
    gsl_vector_free(r);
    gsl_vector_free(xB);
    gsl_matrix_free(Ab_inv);
    gsl_matrix_free(Ab);

    return s;
}