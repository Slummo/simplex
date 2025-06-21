#ifndef PROBLEM_H
#define PROBLEM_H

#define MAX_ROWS 100

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "solution.h"
#include <stdint.h>

typedef struct problem _problem, *problem_t;

// Find problem basis indices with Phase 1 method
uint32_t* problem_find_basis(uint32_t n, uint32_t m, const gsl_matrix* A, const gsl_vector* b, uint32_t* iter_ptr);

// Frees arguments on error
problem_t problem_new(uint32_t n, uint32_t m, uint32_t is_max, gsl_vector* c, gsl_matrix* A, gsl_vector* b,
                      uint32_t* basis);
void problem_print(const problem_t p, const char* name);
void problem_free(problem_t* pp);

solution_t solve(problem_t problem);

#endif