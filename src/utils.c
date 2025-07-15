#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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