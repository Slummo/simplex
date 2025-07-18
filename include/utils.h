#ifndef UTILS_H
#define UTILS_H

#include "utils.h"

#include <stdint.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

int32_t* calculate_nonbasis(int32_t* B, uint32_t n, uint32_t m);

// Print a length‑n vector v in a single row
void print_vector(const char* name, const gsl_vector* v);

// Print an n×m matrix M in a readable form
void print_matrix(const char* name, const gsl_matrix* M);

#endif