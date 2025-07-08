#ifndef SIMPLEX_H
#define SIMPLEX_H

#include "problem.h"
#include "solution.h"
#include "variable.h"
#include <stdint.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

// Find problem basis indices with Phase 1 method
int32_t* simplex_phaseI(problem_t* problem_ptr);

uint32_t simplex_phaseII(uint32_t n, uint32_t m, uint32_t is_max, const gsl_vector* c, const gsl_matrix* A,
                         const gsl_vector* b, int32_t* B, int32_t* N, solution_t* solution_ptr);

#endif