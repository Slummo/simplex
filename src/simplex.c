#include "simplex.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>

// Find problem basis indices with Phase 1 method
int32_t* simplex_phaseI(problem_t* problem_ptr, uint32_t* iter_n_ptr) {
    if (!problem_ptr || !iter_n_ptr) {
        fprintf(stderr, "Some arguments are NULL in simplex_phaseI\n");
        return NULL;
    }

    uint32_t m = problem_m(problem_ptr);
    uint32_t n = problem_n(problem_ptr);
    uint32_t is_max = problem_is_max(problem_ptr);
    gsl_vector* c = problem_c_mut(problem_ptr);
    gsl_matrix* A = problem_A_mut(problem_ptr);
    gsl_vector* b = problem_b_mut(problem_ptr);

    int32_t* B = (int32_t*)malloc(sizeof(int32_t) * n);
    if (!B) {
        fprintf(stderr, "Failed to allocate array for basis indices for PhaseI\n");
        return NULL;
    }

    if (problem_has_feasible_base(problem_ptr, B)) {
        return B;
    }

    int32_t* artificial_B = NULL;
    int32_t* artificial_N = NULL;

    // Augmented capacity for phaseI
    uint32_t variables_num = m + n;
    uint32_t constraints_num = n;

    artificial_B = (int32_t*)malloc(sizeof(int32_t) * n);
    if (!artificial_B) {
        fprintf(stderr, "Failed to allocate array of artificial base indices for PhaseI\n");
        goto fail;
    }

    // Set array
    for (uint32_t i = 0; i < n; i++) {
        artificial_B[i] = (int32_t)(m + i);
    }

    artificial_N = calculate_nonbasis(artificial_B, constraints_num, variables_num);
    if (!artificial_N) {
        fprintf(stderr, "Failed to allocate array of artificial non-base indices for PhaseI\n");
        goto fail;
    }

    // Set the costs for the artificial variables
    for (uint32_t i = m; i < variables_num; i++) {
        gsl_vector_set(c, i, is_max ? -1.0 : 1.0);
    }

    // Add values for artificial variables on each row of the A matrix
    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t j = m; j < variables_num; j++) {
            if (j == m + i) {
                // Artificial variable in column m + i
                gsl_matrix_set(A, i, j, 1.0);
            } else {
                gsl_matrix_set(A, i, j, 0.0);
            }
        }
    }

    var_arr_t* var_arr_ptr = problem_var_arr_mut(problem_ptr);
    variable_t v;
    for (uint32_t i = m; i < variables_num; i++) {
        if (!variable_init_real_positive(&v, 10e9) || !var_arr_push(var_arr_ptr, &v)) {
            goto fail;
        }
    }

    // The PhaseI problem doesn't need another PhaseI, it always
    // has a feasible base and goes straight to PhaseII
    *iter_n_ptr = 0;
    solution_t phaseI_solution;
    if (!simplex_phaseII(constraints_num, variables_num, is_max, c, A, b, artificial_B, artificial_N, &phaseI_solution,
                         iter_n_ptr)) {
        fprintf(stderr, "Failed to run simplex on PhaseI\n");
        goto fail;
    }

    free(artificial_N);
    artificial_N = NULL;

    // No feasible base for original problem
    if (solution_z(&phaseI_solution) < 0) {
        free(B);
        B = NULL;
    } else {
        const gsl_vector* x = solution_x(&phaseI_solution);
        for (uint32_t i = 0; i < n; i++) {
            int32_t var_idx = artificial_B[i];

            // If the variable in basis is artificial (index >= m), make sure its value is 0
            if (var_idx >= (int32_t)m && gsl_vector_get(x, var_idx) > 1e-8) {
                // Infeasible
                free(B);
                B = NULL;
                break;
            }

            B[i] = var_idx;
        }
    }

    free(artificial_B);
    artificial_B = NULL;

    solution_free(&phaseI_solution);

    return B;

fail:
    free(B);
    free(artificial_B);
    free(artificial_N);
    return NULL;
}

uint32_t simplex_phaseII(uint32_t n, uint32_t m, uint32_t is_max, const gsl_vector* c, const gsl_matrix* A,
                         const gsl_vector* b, int32_t* B, int32_t* N, solution_t* solution_ptr, uint32_t* iter_n_ptr) {
    if (!c || !A || !B || !N || !solution_ptr || !iter_n_ptr) {
        fprintf(stderr, "Some arguments are NULL in simplex_phaseII\n");
        return 0;
    }

    gsl_matrix* Ab = gsl_matrix_alloc(n, n);
    gsl_matrix* Ab_inv = NULL;
    gsl_vector* xB = gsl_vector_alloc(n);
    gsl_vector* cb = gsl_vector_alloc(n);
    gsl_vector* cn = gsl_vector_alloc(m - n);
    gsl_vector* r = gsl_vector_alloc(m - n);

    if (!Ab || !xB || !cb || !cn || !r) {
        gsl_matrix_free(Ab);
        gsl_vector_free(xB);
        gsl_vector_free(cb);
        gsl_vector_free(cn);
        gsl_vector_free(r);
        return 0;
    }

    *iter_n_ptr = 0;
    uint32_t unbounded = 0;
    while (1) {
        // Extract Ab matrix and cb vector
        for (uint32_t i = 0; i < n; i++) {
            for (uint32_t j = 0; j < n; j++) {
                // Basis indices must be valid column indices
                gsl_matrix_set(Ab, i, j, gsl_matrix_get(A, i, B[j]));
            }
            double ci = gsl_vector_get(c, B[i]);
            gsl_vector_set(cb, i, is_max ? ci : -ci);
        }

        // Compute the inverse of Ab
        gsl_matrix* old = Ab_inv;
        Ab_inv = inverse(Ab, n);
        gsl_matrix_free(old);

        // Compute xB = Ab_inv * b
        gsl_blas_dgemv(CblasNoTrans, 1.0, Ab_inv, b, 0.0, xB);

        // Compute reduced costs for all nonbasis variables
        for (uint32_t i = 0; i < (m - n); i++) {
            uint32_t j = N[i];

            // cn[i] = cj
            double cj = gsl_vector_get(c, j);
            gsl_vector_set(cn, i, is_max ? cj : -cj);

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

        (*iter_n_ptr)++;
    }

    // Extract optimal solution and value
    if (solution_init(solution_ptr, n, m, unbounded) && !solution_is_unbounded(solution_ptr)) {
        gsl_vector* x = solution_x_mut(solution_ptr);
        double z = 0.0;
        for (uint32_t i = 0; i < n; i++) {
            gsl_vector_set(x, B[i], gsl_vector_get(xB, i));
        }

        gsl_blas_ddot(c, x, &z);
        solution_set_optimal_value(solution_ptr, is_max ? z : -z);
    }

    gsl_vector_free(cb);
    gsl_vector_free(cn);
    gsl_vector_free(r);
    gsl_vector_free(xB);
    gsl_matrix_free(Ab_inv);
    gsl_matrix_free(Ab);

    return 1;
}