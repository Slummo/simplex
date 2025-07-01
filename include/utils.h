#ifndef UTILS_H
#define UTILS_H

#include "rc.h"
#include "variable.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

gsl_vector* vector_read(FILE* stream, char* name, size_t size);
gsl_matrix* matrix_read(FILE* stream, char* name, size_t rows, size_t cols);

gsl_vector* vector_duplicate(const gsl_vector* original);
gsl_matrix* matrix_duplicate(const gsl_matrix* original);

gsl_matrix* inverse(const gsl_matrix* base, size_t size);

#endif