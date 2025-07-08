#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>

uint32_t matrix_init(matrix_t* matrix_ptr, uint32_t row_capacity, uint32_t col_capacity) {
    if (!matrix_ptr) {
        return 0;
    }

    matrix_ptr->data = (double*)calloc(row_capacity * col_capacity, sizeof(double));
    if (!matrix_ptr->data) {
        return 0;
    }

    matrix_ptr->row_capacity = row_capacity;
    matrix_ptr->col_capacity = col_capacity;
    matrix_ptr->stride = col_capacity;

    return 1;
}

uint32_t matrix_from_stream(matrix_t* matrix_ptr, FILE* stream, char* name, uint32_t row_capacity,
                            uint32_t col_capacity, uint32_t rows, uint32_t cols) {
    if (rows > row_capacity || cols > col_capacity) {
        fprintf(stderr, "Requested matrix size exceeds capacity in matrix_from_stream for %s\n", name);
        return 0;
    }

    if (!matrix_init(matrix_ptr, row_capacity, col_capacity)) {
        fprintf(stderr, "Failed to init matrix %s in matrix_from_stream\n", name);
        return 0;
    }

    if (stream == stdin) {
        printf("Enter %ux%u = %u values for matrix %s: ", rows, cols, rows * cols, name);
        fflush(stdout);
    }

    for (uint32_t i = 0; i < rows; i++) {
        for (uint32_t j = 0; j < cols; j++) {
            if (fscanf(stream, "%lf", matrix_ptr->data + i * matrix_ptr->stride + j) != 1) {
                fprintf(stderr, "Failed to read element (%u, %u) of matrix %s\n", i, j, name);
                matrix_free(matrix_ptr);
                return 0;
            }
        }
    }

    return 1;
}

uint32_t matrix_as_gsl_view(matrix_t* matrix_ptr, uint32_t row_offset, uint32_t col_offset, uint32_t rows,
                            uint32_t cols, gsl_matrix_view* view_ptr) {
    if (!matrix_ptr || !matrix_ptr->data) {
        return 0;
    }

    if (row_offset + rows > matrix_ptr->row_capacity || col_offset + cols > matrix_ptr->col_capacity) {
        fprintf(stderr, "View exceeds capacity\n");
        return 0;
    }

    *view_ptr = gsl_matrix_view_array_with_tda(matrix_ptr->data + row_offset * matrix_ptr->stride + col_offset, rows,
                                               cols, matrix_ptr->stride);

    return 1;
}

double matrix_get(const matrix_t* matrix_ptr, uint32_t i, uint32_t j) {
    if (matrix_ptr && i < matrix_ptr->row_capacity && j < matrix_ptr->col_capacity && matrix_ptr->data) {
        return matrix_ptr->data[i * matrix_ptr->stride + j];
    }

    __builtin_trap();
    return 0.0;
}

uint32_t matrix_set(matrix_t* matrix_ptr, uint32_t i, uint32_t j, double value) {
    if (!matrix_ptr || i >= matrix_ptr->row_capacity || j >= matrix_ptr->col_capacity || !matrix_ptr->data) {
        return 0;
    }

    matrix_ptr->data[i * matrix_ptr->stride + j] = value;
    return 1;
}

void matrix_free(matrix_t* matrix_ptr) {
    if (!matrix_ptr) {
        return;
    }

    free(matrix_ptr->data);
    matrix_ptr->data = NULL;
    matrix_ptr->row_capacity = 0;
    matrix_ptr->col_capacity = 0;
    matrix_ptr->stride = 0;
}