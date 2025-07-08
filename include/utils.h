#ifndef UTILS_H
#define UTILS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdint.h>

double* vector_read(FILE* stream, char* name, size_t capacity, size_t size);
double* matrix_read(FILE* stream, char* name, size_t row_capacity, size_t col_capacity, size_t rows, size_t cols);

gsl_vector* vector_duplicate(const gsl_vector* original);
gsl_matrix* matrix_duplicate(const gsl_matrix* original);

gsl_matrix* inverse(const gsl_matrix* base, size_t size);

int32_t* calculate_nonbasis(int32_t* B, uint32_t n, uint32_t m);

#endif