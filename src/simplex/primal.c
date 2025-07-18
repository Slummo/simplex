#include "simplex/primal.h"
#include "simplex/utils.h"
#include "utils.h"

#include <gsl/gsl_linalg.h>

// Find problem basis indices with Phase 1 method
uint32_t simplex_primal_phaseI(uint32_t n, uint32_t m, uint32_t is_max, gsl_vector* c, gsl_matrix* A, gsl_vector* b,
                               int32_t* B, var_arr_t* var_arr_ptr, uint32_t* iter_n_ptr) {
    int32_t* artificial_B = NULL;
    int32_t* artificial_N = NULL;

    // Augmented capacity for phaseI
    uint32_t variables_num = m + n;
    uint32_t constraints_num = n;

    solution_t phaseI_solution = {0};

    uint32_t ret = 1;

    // Create artificial base
    artificial_B = (int32_t*)malloc(sizeof(int32_t) * n);
    if (!artificial_B) {
        fprintf(stderr, "Failed to allocate array of artificial base indices for PhaseI\n");
        goto fail;
    }
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
    for (uint32_t i = 0; i < constraints_num; i++) {
        gsl_matrix_set(A, i, m + i, 1.0);
    }

    variable_t v;
    for (uint32_t i = m; i < variables_num; i++) {
        if (!variable_init_real_positive(&v, 10e9) || !var_arr_push(var_arr_ptr, &v)) {
            goto fail;
        }
    }

    // The PhaseI problem doesn't need another PhaseI, it always
    // has a feasible base and goes straight to PhaseII
    *iter_n_ptr = 0;
    if (!simplex_primal(constraints_num, variables_num, is_max, c, A, b, artificial_B, artificial_N, &phaseI_solution,
                        iter_n_ptr)) {
        fprintf(stderr, "Failed to run primal simplex on PhaseI problem\n");
        goto fail;
    }

    const gsl_vector* x = solution_x(&phaseI_solution);
    for (uint32_t i = 0; i < n; i++) {
        int32_t var_idx = artificial_B[i];
        uint32_t is_var_artificial = var_idx >= (int32_t)m;

        // If the variable in basis is artificial make sure its value is 0
        if (is_var_artificial && gsl_vector_get(x, var_idx) > 1e-8) {
            fprintf(stderr, "No feasible base for original problem\n");
            goto fail;
        }

        B[i] = var_idx;
    }

    goto cleanup;

fail:
    ret = 0;
cleanup:
    free(artificial_B);
    free(artificial_N);
    solution_free(&phaseI_solution);
    return ret;
}

uint32_t simplex_primal(uint32_t n, uint32_t m, uint32_t is_max, const gsl_vector* c, const gsl_matrix* A,
                        const gsl_vector* b, int32_t* B, int32_t* N, solution_t* solution_ptr, uint32_t* iter_n_ptr) {
    if (!c || !A || !B || !N || !solution_ptr || !iter_n_ptr) {
        fprintf(stderr, "Some arguments are NULL in simplex_primal\n");
        return 0;
    }

    uint32_t ret = 1;

    gsl_matrix* AB = gsl_matrix_alloc(n, n);
    gsl_matrix* AB_inv = NULL;
    gsl_vector* xB = gsl_vector_alloc(n);
    gsl_vector* cB = gsl_vector_alloc(n);
    gsl_vector* cN = gsl_vector_alloc(m - n);
    gsl_vector* r = gsl_vector_alloc(m - n);

    if (!AB || !xB || !cB || !cN || !r) {
        goto fail;
    }

    *iter_n_ptr = 0;
    uint32_t unbounded = 0;
    while (1) {
        // Extract AB matrix and cB vector
        extract_basic_objects(n, is_max, B, c, A, cB, AB);

        // Compute AB_inv = AB^-1
        gsl_matrix* old = AB_inv;
        AB_inv = inverse(AB, n);
        gsl_matrix_free(old);
        if (!AB_inv) {
            goto fail;
        }

        // Compute xB = AB_inv * b
        compute_basic_solution(AB_inv, b, xB);

        // For all non-basic variables
        compute_reduced_costs(n, m, is_max, N, c, cB, cN, A, AB_inv, r);

        // Choose the entering variable
        int32_t q = -1;
        for (uint32_t i = 0; i < (m - n); i++) {
            if (gsl_vector_get(r, i) > 0) {
                q = i;
                break;  // Bland's rule: choose first positive reduced cost
            }
        }

        if (q == -1) {
            break;  // Optimal
        }

        // Aq
        gsl_vector* Aq = gsl_vector_alloc(n);
        if (!extract_column(A, (uint32_t)N[q], Aq)) {
            goto fail;
        }

        // Compute direction vector d = -AB_inv * Aq
        gsl_vector* d = gsl_vector_alloc(n);
        gsl_blas_dgemv(CblasNoTrans, -1.0, AB_inv, Aq, 0.0, d);

        // Choose leaving variable
        double min_ratio = 1e20;
        int32_t p = -1;
        for (uint32_t i = 0; i < n; i++) {
            double di = gsl_vector_get(d, i);
            // Only include variables with negative direction coefficient
            if (di < 0.0) {
                double ratio = -gsl_vector_get(xB, i) / di;
                if (ratio < min_ratio) {
                    min_ratio = ratio;
                    p = i;
                }
            }
        }

        // Unbounded
        if (p == -1) {
            unbounded = 1;
            gsl_vector_free(Aq);
            gsl_vector_free(d);
            break;
        }

        pivot(q, p, B, N);

        gsl_vector_free(Aq);
        gsl_vector_free(d);

        (*iter_n_ptr)++;
    }

    // Extract optimal solution and value
    if (solution_init(solution_ptr, n, m + n, unbounded) && !unbounded) {
        extract_optimal(n, is_max, B, xB, c, solution_ptr);
    }

    goto cleanup;

fail:
    ret = 0;

cleanup:
    gsl_matrix_free(AB);
    gsl_matrix_free(AB_inv);
    gsl_vector_free(xB);
    gsl_vector_free(cB);
    gsl_vector_free(cN);
    gsl_vector_free(r);
    return ret;
}