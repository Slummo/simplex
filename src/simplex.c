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

    if (problem_has_primal_feasible_base(problem_ptr, B)) {
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
    for (uint32_t i = 0; i < constraints_num; i++) {
        gsl_matrix_set(A, i, m + i, 1.0);
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
    if (!simplex_primal(constraints_num, variables_num, is_max, c, A, b, artificial_B, artificial_N, &phaseI_solution,
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

void extract_basic_objects(uint32_t n, uint32_t is_max, int32_t* B, const gsl_vector* c, const gsl_matrix* A,
                           gsl_vector* cB, gsl_matrix* AB) {
    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t j = 0; j < n; j++) {
            // Basis indices must be valid column indices
            gsl_matrix_set(AB, i, j, gsl_matrix_get(A, i, B[j]));
        }
        double ci = gsl_vector_get(c, B[i]);
        gsl_vector_set(cB, i, is_max ? ci : -ci);
    }
}

gsl_matrix* inverse(const gsl_matrix* base, size_t size) {
    gsl_matrix* inverse = gsl_matrix_alloc(size, size);
    if (!inverse) {
        return NULL;
    }

    gsl_permutation* p = gsl_permutation_alloc(size);
    if (!p) {
        gsl_matrix_free(inverse);
        return NULL;
    }

    int signum;
    gsl_matrix_memcpy(inverse, base);
    gsl_linalg_LU_decomp(inverse, p, &signum);
    gsl_linalg_LU_invert(inverse, p, inverse);

    gsl_permutation_free(p);

    return inverse;
}

void compute_basic_solution(const gsl_matrix* AB_inv, const gsl_vector* b, gsl_vector* xB) {
    gsl_blas_dgemv(CblasNoTrans, 1.0, AB_inv, b, 0.0, xB);
}

void compute_reduced_costs(uint32_t n, uint32_t m, uint32_t is_max, int32_t* N, const gsl_vector* c, gsl_vector* cB,
                           gsl_vector* cN, const gsl_matrix* A, const gsl_matrix* AB_inv, gsl_vector* r) {
    for (uint32_t i = 0; i < (m - n); i++) {
        uint32_t j = N[i];

        // cN[i] = cj
        double cj = gsl_vector_get(c, j);
        gsl_vector_set(cN, i, is_max ? cj : -cj);

        // Aj
        gsl_vector* Aj = gsl_vector_alloc(n);
        for (uint32_t k = 0; k < n; k++) {
            gsl_vector_set(Aj, k, gsl_matrix_get(A, k, j));
        }

        // AB_inv * Aj
        gsl_vector* Ab_inv_aj = gsl_vector_alloc(n);
        gsl_blas_dgemv(CblasNoTrans, 1.0, AB_inv, Aj, 0.0, Ab_inv_aj);

        // r[i] = rj = cj - cB * AB_inv * Aj
        double res;
        gsl_blas_ddot(cB, Ab_inv_aj, &res);
        gsl_vector_set(r, i, gsl_vector_get(cN, i) - res);

        gsl_vector_free(Aj);
        gsl_vector_free(Ab_inv_aj);
    }
}

uint32_t extract_column(const gsl_matrix* m, uint32_t j, gsl_vector* col) {
    if (!m || !col) {
        return 0;
    }

    for (uint32_t i = 0; i < col->size; i++) {
        gsl_vector_set(col, i, gsl_matrix_get(m, i, j));
    }

    return 1;
}

uint32_t extract_row(const gsl_matrix* m, uint32_t i, gsl_vector* row) {
    if (!m || !row) {
        return 0;
    }

    for (uint32_t j = 0; j < row->size; j++) {
        gsl_vector_set(row, j, gsl_matrix_get(m, i, j));
    }

    return 1;
}

void pivot(int32_t entering, int32_t leaving, int32_t* B, int32_t* N) {
    uint32_t tmp = B[leaving];
    B[leaving] = N[entering];
    N[entering] = tmp;
}

void extract_optimal(uint32_t n, uint32_t is_max, int32_t* B, gsl_vector* xB, const gsl_vector* c,
                     solution_t* solution_ptr) {
    gsl_vector* x = solution_x_mut(solution_ptr);
    double z = 0.0;
    for (uint32_t i = 0; i < n; i++) {
        gsl_vector_set(x, B[i], gsl_vector_get(xB, i));
    }

    gsl_blas_ddot(c, x, &z);
    solution_set_optimal_value(solution_ptr, is_max ? z : -z);
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
        gsl_matrix_scale(AB_inv, -1.0);
        gsl_blas_dgemv(CblasNoTrans, 1.0, AB_inv, Aq, 0.0, d);

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