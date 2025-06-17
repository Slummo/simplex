#ifndef UTILS_H
#define UTILS_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

gsl_vector* read_vector(FILE* stream, char* name, size_t size);
gsl_matrix* read_matrix(FILE* stream, char* name, size_t rows, size_t cols);

gsl_matrix* inverse(const gsl_matrix* base, size_t size);

#endif