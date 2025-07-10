#ifndef UTILS_H
#define UTILS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdint.h>

gsl_vector* gsl_vector_from_stream(FILE* stream, char* name, uint32_t capacity, uint32_t size);
gsl_matrix* gsl_matrix_from_stream(FILE* stream, char* name, uint32_t row_capacity, uint32_t col_capacity,
                                   uint32_t rows, uint32_t cols);

gsl_vector* vector_duplicate(const gsl_vector* original);
gsl_matrix* matrix_duplicate(const gsl_matrix* original);

gsl_matrix* inverse(const gsl_matrix* base, size_t size);

int32_t* calculate_nonbasis(int32_t* B, uint32_t n, uint32_t m);

#endif