#include "simplex/utils.h"

#include <gsl/gsl_linalg.h>

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
