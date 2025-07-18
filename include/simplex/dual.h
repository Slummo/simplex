#ifndef DUAL_H
#define DUAL_H

#include "solution.h"

#include <gsl/gsl_matrix.h>

uint32_t simplex_dual(uint32_t n, uint32_t m, uint32_t is_max, const gsl_vector* c, const gsl_matrix* A,
                      const gsl_vector* b, int32_t* B, int32_t* N, solution_t* solution_ptr, uint32_t* iter_n_ptr);

#endif