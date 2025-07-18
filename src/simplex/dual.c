#include "simplex/dual.h"
#include "simplex/utils.h"

#include <gsl/gsl_linalg.h>

uint32_t simplex_dual(uint32_t n, uint32_t m, uint32_t is_max, const gsl_vector* c, const gsl_matrix* A,
                      const gsl_vector* b, int32_t* B, int32_t* N, solution_t* solution_ptr, uint32_t* iter_n_ptr) {
    if (!c || !A || !B || !N || !solution_ptr || !iter_n_ptr) {
        fprintf(stderr, "Some arguments are NULL in simplex_dual\n");
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

        // Choose leaving basic variable (primal-infeasible)
        int32_t p = -1;
        double most_negative = -1e-8;
        for (uint32_t i = 0; i < n; i++) {
            double xi = gsl_vector_get(xB, i);
            if (xi < most_negative) {
                most_negative = xi;
                p = i;
            }
        }

        if (p == -1) {
            break;  // Primal feasible, so optimal
        }

        // Extract leaving row from AB_inv
        gsl_vector* AB_inv_p = gsl_vector_alloc(n);
        if (!extract_row(AB_inv, (uint32_t)p, AB_inv_p)) {
            goto fail;
        }

        // Choose entering variable by computing alpha_pj = AB_inv_p * Aj for each non-basic variable
        double min_ratio = 1e20;
        int32_t q = -1;
        gsl_vector* Aj = gsl_vector_alloc(n);
        for (uint32_t i = 0; i < m - n; i++) {
            uint32_t j = N[i];

            // Aj
            if (!extract_column(A, j, Aj)) {
                gsl_vector_free(AB_inv_p);
                gsl_vector_free(Aj);
                goto fail;
            }

            double alpha_pj;
            gsl_blas_ddot(AB_inv_p, Aj, &alpha_pj);

            // Only include positive ones
            if (alpha_pj > 0.0) {
                double ratio = -gsl_vector_get(r, j) / alpha_pj;
                if (ratio < min_ratio) {
                    min_ratio = ratio;
                    q = i;
                }
            }
        }
        gsl_vector_free(Aj);
        gsl_vector_free(AB_inv_p);

        if (q == -1) {
            unbounded = 1;
            break;
        }

        pivot(q, p, B, N);

        (*iter_n_ptr)++;
    }

    // Extract optimal solution and value
    if (solution_init(solution_ptr, n, m, unbounded) && !unbounded) {
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