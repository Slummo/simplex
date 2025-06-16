#include "problem.h"
#include <gsl/gsl_linalg.h>
#include <stdlib.h>
#include <string.h>

struct problem {
    int n;          // Number of constraints
    int m;          // Number of variables
    gsl_vector* c;  // Reduced costs (m)
    gsl_matrix* A;  // Constraints matrix (n x m)
    gsl_vector* b;  // RHS (n)
    int* basis;     // Indices of basic variables (size n)
    int* nonbasis;  // Indices of nonbasic variables (size m-n)
};

problem_t problem_new(int n, int m, gsl_vector* c, gsl_matrix* A, gsl_vector* b, int* basis) {
    problem_t p = (problem_t)malloc(sizeof(_problem));
    if (!p) {
        return NULL;
    }

    p->n = n;
    p->m = m;
    p->c = c;
    p->A = A;
    p->b = b;
    p->basis = basis;

    p->nonbasis = (int*)malloc(sizeof(int) * (m - n));
    if (!p->nonbasis) {
        return NULL;
    }

    // Fill non-basis from remaining indices
    int* used = calloc(m, sizeof(int));
    for (int i = 0; i < n; i++) {
        if (used[basis[i]]) {
            fprintf(stderr, "Duplicate basis index %d\n", basis[i]);
            problem_free(&p);
            return NULL;
        }
        used[basis[i]] = 1;
    }

    // Fill non-basis with those not used
    int ni = 0;
    for (int i = 0; i < m; i++) {
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

struct solution {
    int n;             // Number of variables
    gsl_vector* x;     // Optimal solution (n)
    double z;          // Optimal value
    int is_unbounded;  // Boolean value to know if unbounded
    int n_iter;        // Number of iterations needed to solve
};

solution_t solution_new(int n, gsl_vector* x, double z, int is_unbounded, int n_iter) {
    solution_t s = (solution_t)malloc(sizeof(_solution));
    if (!s) {
        return NULL;
    }

    s->n = n;
    s->x = x;
    s->z = z;
    s->is_unbounded = is_unbounded;
    s->n_iter = n_iter;

    return s;
}

void solution_print(solution_t s) {
    if (!s) {
        return;
    }

    printf("\n==================SOLUTION==================\n");
    if (s->is_unbounded) {
        printf("infinite\n");
    } else {
        printf("Optimal found in %d iterations\nz*: %lf\nx*: (", s->n_iter, s->z);
    }

    if (!s->is_unbounded) {
        for (int i = 0; i < s->n; i++) {
            printf("%.3lf", gsl_vector_get(s->x, i));
            if (i < s->n - 1) {
                printf(", ");
            }
        }
        puts(")");
    }
}

void solution_free(solution_t* sp) {
    if (!sp || !*sp) {
        return;
    }

    gsl_vector_free((*sp)->x);
    free(*sp);
    *sp = NULL;
}

solution_t solve(problem_t problem) {
    if (!problem) {
        return NULL;
    }

    int n = problem->n;
    int m = problem->m;
    gsl_matrix* A = problem->A;
    gsl_vector* b = problem->b;
    gsl_vector* c = problem->c;
    int* B = problem->basis;
    int* N = problem->nonbasis;

    gsl_matrix* Ab = gsl_matrix_alloc(n, n);
    gsl_matrix* Ab_inv = gsl_matrix_alloc(n, n);
    gsl_vector* xB = gsl_vector_alloc(n);
    gsl_vector* cb = gsl_vector_alloc(n);
    gsl_vector* cn = gsl_vector_alloc(m - n);
    gsl_vector* r = gsl_vector_alloc(m - n);

    int iter = 0;
    int unbounded = 0;
    while (1) {
        // Extract Ab matrix
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                gsl_matrix_set(Ab, i, j, gsl_matrix_get(A, i, B[j]));
            }
            gsl_vector_set(cb, i, gsl_vector_get(c, B[i]));
        }

        // Compute the inverse of Ab
        gsl_permutation* p = gsl_permutation_alloc(n);
        int signum;
        gsl_matrix_memcpy(Ab_inv, Ab);
        gsl_linalg_LU_decomp(Ab_inv, p, &signum);
        gsl_linalg_LU_invert(Ab_inv, p, Ab_inv);
        gsl_permutation_free(p);

        // Compute xB = Ab_inv * b
        gsl_blas_dgemv(CblasNoTrans, 1.0, Ab_inv, b, 0.0, xB);

        // Compute reduced costs for all nonbasis variables
        for (int i = 0; i < m - n; i++) {
            int j = N[i];

            // cn[i] = cj
            gsl_vector_set(cn, i, gsl_vector_get(c, j));

            // aj
            gsl_vector* aj = gsl_vector_alloc(n);
            for (int k = 0; k < n; k++) {
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
        for (int i = 0; i < m - n; i++) {
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
        for (int k = 0; k < n; k++) {
            gsl_vector_set(aj, k, gsl_matrix_get(A, k, j));
        }

        // Compute direction vector d = -Ab_inv * aj
        gsl_vector* d = gsl_vector_alloc(n);
        gsl_matrix_scale(Ab_inv, -1.0);
        gsl_blas_dgemv(CblasNoTrans, 1.0, Ab_inv, aj, 0.0, d);

        // Choose leaving variable
        double min_ratio = 1e20;
        int leaving = -1;
        for (int i = 0; i < n; i++) {
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

    // Extract optimal solution and value
    if (!s->is_unbounded) {
        for (int i = 0; i < n; i++) {
            gsl_vector_set(s->x, B[i], gsl_vector_get(xB, i));
        }
        gsl_blas_ddot(c, s->x, &s->z);
    }

    gsl_vector_free(cb);
    gsl_vector_free(cn);
    gsl_vector_free(r);
    gsl_vector_free(xB);
    gsl_matrix_free(Ab_inv);
    gsl_matrix_free(Ab);

    return s;
}