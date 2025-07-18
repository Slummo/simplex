#ifndef PRIMAL_H
#define PRIMAL_H

#include "problem.h"

// Find problem basis indices with Phase 1 method
uint32_t simplex_primal_phaseI(uint32_t n, uint32_t m, uint32_t is_max, gsl_vector* c, gsl_matrix* A, gsl_vector* b,
                               int32_t* B, var_arr_t* var_arr_ptr, uint32_t* iter_n_ptr);

uint32_t simplex_primal(uint32_t n, uint32_t m, uint32_t is_max, const gsl_vector* c, const gsl_matrix* A,
                        const gsl_vector* b, int32_t* B, int32_t* N, solution_t* solution_ptr, uint32_t* iter_n_ptr);

#endif