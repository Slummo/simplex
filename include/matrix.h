#ifndef MATRIX_H
#define MATRIX_H

#include <stdint.h>
#include <gsl/gsl_matrix.h>

typedef struct matrix {
    double* data;
    uint32_t row_capacity;
    uint32_t col_capacity;
    uint32_t stride;  // Number of elements between rows (columns)
} matrix_t;

uint32_t matrix_init(matrix_t* matrix_ptr, uint32_t row_capacity, uint32_t col_capacity);
uint32_t matrix_from_stream(matrix_t* matrix_ptr, FILE* stream, char* name, uint32_t row_capacity,
                            uint32_t col_capacity, uint32_t rows, uint32_t cols);
uint32_t matrix_as_gsl_view(matrix_t* matrix_ptr, uint32_t row_offset, uint32_t col_offset, uint32_t rows,
                            uint32_t cols, gsl_matrix_view* view_ptr);
double matrix_get(const matrix_t* matrix_ptr, uint32_t i, uint32_t j);
uint32_t matrix_set(matrix_t* matrix_ptr, uint32_t i, uint32_t j, double value);
void matrix_free(matrix_t* matrix_ptr);

#endif