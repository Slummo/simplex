#include "utils.h"
#include <gsl/gsl_linalg.h>

gsl_vector* vector_read(FILE* stream, char* name, size_t size) {
    gsl_vector* v = gsl_vector_alloc(size);
    if (!v) {
        fprintf(stderr, "Failed to allocate vector %s in vector_read\n", name);
        return NULL;
    }

    if (stream == stdin) {
        printf("Enter %zu values for vector %s: ", size, name);
        fflush(stdout);
    }

    if (gsl_vector_fscanf(stream, v) != 0) {
        fprintf(stderr, "Failed to read vector %s in vector_read\n", name);
        gsl_vector_free(v);
        return NULL;
    }

    return v;
}

gsl_matrix* matrix_read(FILE* stream, char* name, size_t rows, size_t cols) {
    gsl_matrix* m = gsl_matrix_alloc(rows, cols);
    if (!m) {
        fprintf(stderr, "Failed to allocate matrix %s in matrix_read\n", name);
        return NULL;
    }

    if (stream == stdin) {
        printf("Enter %zux%zu = %zu values for matrix %s: ", rows, cols, rows * cols, name);
        fflush(stdout);
    }

    if (gsl_matrix_fscanf(stream, m) != 0) {
        fprintf(stderr, "Failed to read matrix %s in matrix_read\n", name);
        gsl_matrix_free(m);
        return NULL;
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