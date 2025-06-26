#include "simplex.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_linalg.h>

// Find problem basis indices with Phase 1 method
int32_t* simplex_phaseI(uint32_t n, uint32_t m, const gsl_matrix* A, const gsl_vector* b, const variable_t* variables,
                        uint32_t* pI_iter_ptr) {
    if (!A || !b || !pI_iter_ptr) {
        return NULL;
    }

    int32_t* basis = NULL;
    gsl_vector* c = NULL;
    int32_t* artificial_indices = NULL;
    gsl_matrix* A2 = NULL;
    gsl_vector* b2 = NULL;
    variable_t* variables2 = NULL;
    problem_t phaseI = NULL;

    // Try first to identify an identity matrix
    basis = (int32_t*)malloc(sizeof(int32_t) * n);
    if (!basis) {
        fprintf(stderr, "Failed to allocate array for basis indices for PhaseI\n");
        goto fail;
    }

    // Initialize basis indices to -1
    for (uint32_t i = 0; i < n; i++) {
        basis[i] = -1;
    }

    // Check for identity
    for (uint32_t j = 0; j < m; j++) {
        int32_t pivot_row = -1;
        uint32_t valid = 1;
        for (uint32_t i = 0; i < n; i++) {
            double aij = gsl_matrix_get(A, i, j);
            if (aij == (double)1) {
                if (pivot_row == -1) {
                    pivot_row = (int32_t)i;
                } else {
                    valid = 0;
                    break;
                }
            } else if (aij != (double)0) {
                valid = 0;
                break;
            }
        }
        if (valid && pivot_row != -1 && basis[pivot_row] == -1) {
            basis[pivot_row] = j;
        }
    }

    // Check if base is feasible
    uint32_t feasible = 1;
    for (uint32_t i = 0; i < n; i++) {
        if (basis[i] == -1 || gsl_vector_get(b, i) < 0) {
            feasible = 0;
            break;
        }
    }

    if (feasible) {
        pI_iter_ptr = 0;
        return basis;
    }

    // If not, try with Phase I method.
    memset(basis, 0, sizeof(int32_t) * n);

    // One artificial variable per constraint
    c = gsl_vector_alloc(m + n);
    if (!c) {
        fprintf(stderr, "Failed to alloc c for PhaseI\n");
        goto fail;
    }

    // Set vector c
    for (uint32_t i = 0; i < m + n; i++) {
        gsl_vector_set(c, i, i < m ? 0.0 : -1.0);
    }

    artificial_indices = (int32_t*)malloc(sizeof(int32_t) * n);
    if (!artificial_indices) {
        fprintf(stderr, "Failed to allocate array for artificial indices for PhaseI\n");
        goto fail;
    }

    // Set array
    for (uint32_t i = 0; i < n; i++) {
        artificial_indices[i] = (int32_t)(m + i);
    }

    // Augmented matrix with m + n columns (one column for each artificial variable)
    A2 = gsl_matrix_alloc(n, m + n);
    if (!A2) {
        fprintf(stderr, "Failed to alloc A2 for PhaseI\n");
        goto fail;
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

    // Duplicate b and variables
    b2 = vector_duplicate(b);
    if (!b2) {
        fprintf(stderr, "Failed to duplicate b for PhaseI\n");
        goto fail;
    }

    variables2 = (variable_t*)malloc(sizeof(variable_t) * m);
    if (!variables2) {
        fprintf(stderr, "Failed to duplicate variables for PhaseI\n");
        goto fail;
    }

    for (uint32_t i = 0; i < m; i++) {
        variables2[i] = variables[i];
    }

    // The PhaseI problem doesn't need another PhaseI, it always
    // has a feasible base and goes straight to PhaseII
    phaseI = problem_new(n, m + n, 1, c, A2, b2, artificial_indices, 0, variables2);
    if (!phaseI) {
        fprintf(stderr, "Failed to create PhaseI problem\n");
        goto fail;
    }

    solution_t pII_s = simplex_phaseII(phaseI);
    if (!pII_s) {
        fprintf(stderr, "Failed to create solution for PhaseI problem\n");
        goto fail;
    }

    // No feasible base for original problem
    if (solution_z(pII_s) < 0) {
        free(basis);
        basis = NULL;
    } else {
        const gsl_vector* x = solution_x(pII_s);
        const int32_t* final_basis = solution_basis(pII_s);
        for (uint32_t i = 0; i < n; i++) {
            // If the variable in basis is artificial (index >= m), make sure its value is 0
            if (final_basis[i] >= (int32_t)m && gsl_vector_get(x, final_basis[i]) > 0) {
                // Infeasible
                free(basis);
                basis = NULL;
                break;
            }

            basis[i] = final_basis[i];
        }
    }

    *pI_iter_ptr = solution_pII_iterations(pII_s);

    gsl_vector_free(c);
    free(artificial_indices);
    gsl_matrix_free(A2);
    gsl_vector_free(b2);
    variables_arr_free(&variables2, m);
    problem_free(&phaseI);
    solution_free(&pII_s);

    return basis;

fail:
    free(basis);
    gsl_vector_free(c);
    free(artificial_indices);
    gsl_matrix_free(A2);
    gsl_vector_free(b2);
    variables_arr_free(&variables2, m);
    problem_free(&phaseI);
    return NULL;
}

// simplex_phaseII method on linear problem p
solution_t simplex_phaseII(const problem_t p) {
    if (!p) {
        fprintf(stderr, "problem is NULL in simplex_phaseII\n");
        return NULL;
    }

    uint32_t n = problem_n(p);
    uint32_t m = problem_m(p);
    const gsl_matrix* A = problem_A(p);
    const gsl_vector* b = problem_b(p);
    gsl_vector* c = vector_duplicate(problem_c(p));
    if (!c) {
        return NULL;
    }
    if (!problem_is_max(p)) {
        gsl_vector_scale(c, -1.0);
    }

    int32_t* basis = problem_basis_mut(p);
    int32_t* nonbasis = problem_nonbasis_mut(p);

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
        return NULL;
    }

    uint32_t pII_iter = 0;
    uint32_t unbounded = 0;
    while (1) {
        // Extract Ab matrix
        for (uint32_t i = 0; i < n; i++) {
            for (uint32_t j = 0; j < n; j++) {
                // basis indices must be valid column indices
                gsl_matrix_set(Ab, i, j, gsl_matrix_get(A, i, basis[j]));
            }
            gsl_vector_set(cb, i, gsl_vector_get(c, basis[i]));
        }

        // Compute the inverse of Ab
        Ab_inv = inverse(Ab, n);

        // Compute xB = Ab_inv * b
        gsl_blas_dgemv(CblasNoTrans, 1.0, Ab_inv, b, 0.0, xB);

        // Compute reduced costs for all nonbasis variables
        for (uint32_t i = 0; i < (m - n); i++) {
            uint32_t j = nonbasis[i];

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
        uint32_t j = nonbasis[entering];
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
        uint32_t tmp = basis[leaving];
        basis[leaving] = nonbasis[entering];
        nonbasis[entering] = tmp;

        gsl_vector_free(aj);
        gsl_vector_free(d);
        pII_iter++;
    }

    solution_t s = solution_new(m, gsl_vector_calloc(m), basis, unbounded, problem_pI_iter(p), pII_iter);

    // Extract optimal solution and value
    if (s && !solution_is_unbounded(s)) {
        gsl_vector* v = solution_x_mut(s);
        double z = 0.0;
        for (uint32_t i = 0; i < n; i++) {
            gsl_vector_set(v, basis[i], gsl_vector_get(xB, i));
        }
        gsl_blas_ddot(c, v, &z);
        solution_set_optimal_value(s, z);
    }

    gsl_vector_free(c);
    gsl_vector_free(cb);
    gsl_vector_free(cn);
    gsl_vector_free(r);
    gsl_vector_free(xB);
    gsl_matrix_free(Ab_inv);
    gsl_matrix_free(Ab);

    return s;
}