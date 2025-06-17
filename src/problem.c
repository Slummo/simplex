#include "problem.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_linalg.h>

struct problem {
    unsigned int n;          // Number of constraints
    unsigned int m;          // Number of variables
    int is_max;              // Boolean value to know if its a maximization problem
    gsl_vector* c;           // Reduced costs (m)
    gsl_matrix* A;           // Constraints matrix (n x m)
    gsl_vector* b;           // RHS (n)
    unsigned int* basis;     // Indices of basic variables (size n)
    unsigned int* nonbasis;  // Indices of nonbasic variables (size m-n)
};

problem_t problem_new(size_t n, size_t m, int is_max, gsl_vector* c, gsl_matrix* A, gsl_vector* b,
                      unsigned int* basis) {
    problem_t p = (problem_t)malloc(sizeof(_problem));
    if (!p) {
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

    p->nonbasis = (unsigned int*)malloc(sizeof(unsigned int) * (m - n));
    if (!p->nonbasis) {
        return NULL;
    }

    // Fill non-basis from remaining indices
    unsigned int* used = calloc(m, sizeof(unsigned int));
    for (int i = 0; i < (int)n; i++) {
        if (used[basis[i]]) {
            fprintf(stderr, "Duplicate basis index %d\n", basis[i]);
            problem_free(&p);
            return NULL;
        }
        used[basis[i]] = 1;
    }

    // Fill non-basis with those not used
    int ni = 0;
    for (int i = 0; i < (int)m; i++) {
        if (!used[i]) {
            p->nonbasis[ni++] = i;
        }
    }
    free(used);

    return p;
}

void problem_free(problem_t* pp) {
    if (!pp || !*pp) {
        return;
    }

    gsl_vector_free((*pp)->c);
    gsl_matrix_free((*pp)->A);
    gsl_vector_free((*pp)->b);
    free((*pp)->basis);
    free((*pp)->nonbasis);
    free(*pp);
    *pp = NULL;
}

solution_t solve(problem_t problem) {
    if (!problem) {
        return NULL;
    }

    size_t n = problem->n;
    size_t m = problem->m;
    gsl_matrix* A = problem->A;
    gsl_vector* b = problem->b;
    gsl_vector* c = problem->c;
    unsigned int* B = problem->basis;
    unsigned int* N = problem->nonbasis;

    gsl_matrix* Ab = gsl_matrix_alloc(n, n);
    gsl_matrix* Ab_inv = NULL;
    gsl_vector* xB = gsl_vector_alloc(n);
    gsl_vector* cb = gsl_vector_alloc(n);
    gsl_vector* cn = gsl_vector_alloc(m - n);
    gsl_vector* r = gsl_vector_alloc(m - n);

    int iter = 0;
    int unbounded = 0;
    while (1) {
        // Extract Ab matrix
        for (int i = 0; i < (int)n; i++) {
            for (int j = 0; j < (int)n; j++) {
                gsl_matrix_set(Ab, i, j, gsl_matrix_get(A, i, B[j]));
            }
            gsl_vector_set(cb, i, gsl_vector_get(c, B[i]));
        }

        // Compute the inverse of Ab
        Ab_inv = inverse(Ab, n);

        // Compute xB = Ab_inv * b
        gsl_blas_dgemv(CblasNoTrans, 1.0, Ab_inv, b, 0.0, xB);

        // Compute reduced costs for all nonbasis variables
        for (int i = 0; i < (int)(m - n); i++) {
            int j = N[i];

            // cn[i] = cj
            gsl_vector_set(cn, i, gsl_vector_get(c, j));

            // aj
            gsl_vector* aj = gsl_vector_alloc(n);
            for (int k = 0; k < (int)n; k++) {
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
        int entering = -1;
        for (int i = 0; i < (int)(m - n); i++) {
            if (gsl_vector_get(r, i) > 0) {
                entering = i;
                break;  // Bland's rule: choose first positive reduced cost
            }
        }

        if (entering == -1) {
            break;  // Optimal
        }

        // Extract column j of A
        int j = N[entering];
        gsl_vector* aj = gsl_vector_alloc(n);
        for (int k = 0; k < (int)n; k++) {
            gsl_vector_set(aj, k, gsl_matrix_get(A, k, j));
        }

        // Compute direction vector d = -Ab_inv * aj
        gsl_vector* d = gsl_vector_alloc(n);
        gsl_matrix_scale(Ab_inv, -1.0);
        gsl_blas_dgemv(CblasNoTrans, 1.0, Ab_inv, aj, 0.0, d);

        // Choose leaving variable
        double min_ratio = 1e20;
        int leaving = -1;
        for (int i = 0; i < (int)n; i++) {
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
        int tmp = B[leaving];
        B[leaving] = N[entering];
        N[entering] = tmp;

        gsl_vector_free(aj);
        gsl_vector_free(d);
        iter++;
    }

    solution_t s = solution_new(m, gsl_vector_calloc(m), 0.0, unbounded, iter);
    if (!s) {
        gsl_vector_free(cb);
        gsl_vector_free(cn);
        gsl_vector_free(r);
        gsl_vector_free(xB);
        gsl_matrix_free(Ab_inv);
        gsl_matrix_free(Ab);
        return NULL;
    }

    // Extract optimal solution and value
    if (!is_solution_unbounded(s)) {
        gsl_vector* v = solution_optimal_vec_ptr(s);
        for (int i = 0; i < (int)n; i++) {
            gsl_vector_set(v, B[i], gsl_vector_get(xB, i));
        }
        gsl_blas_ddot(c, v, solution_optimal_value_ptr(s));
    }

    gsl_vector_free(cb);
    gsl_vector_free(cn);
    gsl_vector_free(r);
    gsl_vector_free(xB);
    gsl_matrix_free(Ab_inv);
    gsl_matrix_free(Ab);

    return s;
}