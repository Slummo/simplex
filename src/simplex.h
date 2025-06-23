#ifndef SIMPLEX_H
#define SIMPLEX_H

#include "problem.h"
#include "solution.h"
#include <stdint.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

// Find problem basis indices with Phase 1 method
int32_t* simplex_phaseI(uint32_t n, uint32_t m, const gsl_matrix* A, const gsl_vector* b, uint32_t* pI_iter_ptr);

// simplex_phaseII method on linear problem p
solution_t simplex_phaseII(const problem_t problem, uint32_t pII_iter);

#endif