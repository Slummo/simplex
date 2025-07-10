#include "utils.h"
#include <gsl/gsl_linalg.h>
#include <string.h>

gsl_vector* gsl_vector_from_stream(FILE* stream, char* name, uint32_t capacity, uint32_t size) {
    if (size > capacity) {
        fprintf(stderr, "Requested size exceeds capacity in gsl_vector_from_stream for %s\n", name);
        return NULL;
    }

    gsl_vector* v = gsl_vector_calloc(capacity);
    if (!v) {
        fprintf(stderr, "Failed to allocate gsl_vector %s in gsl_vector_from_stream\n", name);
        return NULL;
    }

    if (stream == stdin) {
        printf("Enter %u values for vector %s: ", size, name);
        fflush(stdout);
    }

    for (uint32_t i = 0; i < size; i++) {
        if (fscanf(stream, "%lf", v->data + i) != 1) {
            fprintf(stderr, "Failed to read element %u of gsl_vector %s\n", i, name);
            gsl_vector_free(v);
            return NULL;
        }
    }

    return v;
}

gsl_matrix* gsl_matrix_from_stream(FILE* stream, char* name, uint32_t row_capacity, uint32_t col_capacity,
                                   uint32_t rows, uint32_t cols) {
    if (rows > row_capacity || cols > col_capacity) {
        fprintf(stderr, "Requested size exceeds capacity in gsl_matrix_from_stream for %s\n", name);
        return NULL;
    }

    gsl_matrix* m = gsl_matrix_calloc(row_capacity, col_capacity);
    if (!m) {
        fprintf(stderr, "Failed to allocate gsl_matrix %s in gsl_matrix_from_stream\n", name);
        return NULL;
    }

    if (stream == stdin) {
        printf("Enter %ux%u = %u values for matrix %s: ", rows, cols, rows * cols, name);
        fflush(stdout);
    }

    for (uint32_t i = 0; i < rows; i++) {
        for (uint32_t j = 0; j < cols; j++) {
            if (fscanf(stream, "%lf", m->data + i * m->tda + j) != 1) {
                fprintf(stderr, "Failed to read element (%u, %u) of gsl_matrix %s\n", i, j, name);
                gsl_matrix_free(m);
                return NULL;
            }
        }
    }

    return m;
}

gsl_vector* vector_duplicate(const gsl_vector* original) {
    if (!original) {
        fprintf(stderr, "Original vector is NULL in vector_duplicate\n");
        return NULL;
    }

    gsl_vector* new = gsl_vector_alloc(original->size);
    if (!new) {
        fprintf(stderr, "Failed to allocate vector in vector_duplicate\n");
        return NULL;
    }

    gsl_vector_memcpy(new, original);

    return new;
}

gsl_matrix* matrix_duplicate(const gsl_matrix* original) {
    if (!original) {
        fprintf(stderr, "Original matrix is NULL in matrix_duplicate\n");
        return NULL;
    }

    gsl_matrix* new = gsl_matrix_alloc(original->size1, original->size2);
    if (!new) {
        fprintf(stderr, "Failed to allocate matrix in matrix_duplicate\n");
        return NULL;
    }

    gsl_matrix_memcpy(new, original);

    return new;
}

gsl_matrix* inverse(const gsl_matrix* base, size_t size) {
    gsl_matrix* inverse = gsl_matrix_alloc(size, size);
    if (!inverse) {
        return NULL;
    }

    gsl_permutation* p = gsl_permutation_alloc(size);
    if (!p) {
        gsl_matrix_free(inverse);
    }

    int signum;
    gsl_matrix_memcpy(inverse, base);
    gsl_linalg_LU_decomp(inverse, p, &signum);
    gsl_linalg_LU_invert(inverse, p, inverse);

    gsl_permutation_free(p);

    return inverse;
}

int32_t* calculate_nonbasis(int32_t* B, uint32_t n, uint32_t m) {
    if (!B) {
        return NULL;
    }

    // Create nonbasis indices array
    int32_t* used = (int32_t*)malloc(sizeof(int32_t) * m);
    if (!used) {
        return NULL;
    }
    memset(used, 0, sizeof(int32_t) * m);

    for (uint32_t i = 0; i < n; i++) {
        if (B[i] >= (int32_t)m || used[B[i]]) {
            fprintf(stderr, "Duplicate or invalid basis index %u\n", B[i]);
            free(used);
            return NULL;
        }
        used[B[i]] = 1;
    }

    int32_t* N = (int32_t*)malloc(sizeof(int32_t) * (m - n));
    if (!N) {
        free(used);
        return NULL;
    }

    uint32_t nonbasis_count = 0;
    for (uint32_t i = 0; i < m; i++) {
        if (!used[i]) {
            N[nonbasis_count++] = i;
        }
    }

    free(used);
    return N;
}