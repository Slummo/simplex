#include "utils.h"

gsl_vector* read_vector(FILE* stream, char* name, size_t size) {
    gsl_vector* v = gsl_vector_alloc(size);
    if (!v) {
        fprintf(stderr, "Failed to allocate vector %s\n", name);
        return NULL;
    }

    printf("Enter %zu values for vector %s: ", size, name);
    if (gsl_vector_fscanf(stream, v) != 0) {
        fprintf(stderr, "Failed to read vector %s\n", name);
        gsl_vector_free(v);
        return NULL;
    }

    return v;
}

gsl_matrix* read_matrix(FILE* stream, char* name, size_t rows, size_t cols) {
    gsl_matrix* m = gsl_matrix_alloc(rows, cols);
    if (!m) {
        fprintf(stderr, "Failed to allocate matrix %s\n", name);
        return NULL;
    }

    printf("Enter %zux%zu = %zu values for matrix %s: ", rows, cols, rows * cols, name);
    if (gsl_matrix_fscanf(stream, m) != 0) {
        fprintf(stderr, "Failed to read matrix %s\n", name);
        gsl_matrix_free(m);
        return NULL;
    }

    return m;
}

gsl_matrix* inverse(const gsl_matrix* base, size_t size) {
    gsl_matrix* inverse = gsl_matrix_alloc(size, size);
    gsl_permutation* p = gsl_permutation_alloc(size);
    if (!inverse || !p) {
        return NULL;
    }

    int signum;
    gsl_matrix_memcpy(inverse, base);
    gsl_linalg_LU_decomp(inverse, p, &signum);
    gsl_linalg_LU_invert(inverse, p, inverse);
    gsl_permutation_free(p);

    return inverse;
}