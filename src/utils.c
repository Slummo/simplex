#include "utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int32_t* calculate_nonbasis(int32_t* B, uint32_t n, uint32_t m) {
    if (!B) {
        return NULL;
    }

    // Create nonbasis indices array
    int32_t* used = (int32_t*)calloc(m, sizeof(int32_t));
    if (!used) {
        return NULL;
    }

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

// Print a length‑n vector v in a single row.
void print_vector(const char* name, const gsl_vector* v) {
    size_t n = v->size;
    printf("%s = [", name);
    for (size_t i = 0; i < n; i++) {
        printf(" %.6g", gsl_vector_get(v, i));
        if (i + 1 < n)
            putchar(',');
    }
    printf(" ]\n");
}

// Print an n×m matrix M in a readable form.
void print_matrix(const char* name, const gsl_matrix* M) {
    size_t n = M->size1;
    size_t m = M->size2;
    printf("%s = [\n", name);
    for (size_t i = 0; i < n; i++) {
        printf("  [");
        for (size_t j = 0; j < m; j++) {
            printf(" % .6g", gsl_matrix_get(M, i, j));
            if (j + 1 < m)
                putchar(',');
        }
        printf(" ]\n");
    }
    printf("]\n");
}
