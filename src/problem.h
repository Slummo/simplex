#ifndef PROBLEM_H
#define PROBLEM_H

#define MAX_ROWS 100

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "solution.h"

typedef struct problem _problem, *problem_t;
problem_t problem_new(size_t n, size_t m, int is_max, gsl_vector* c, gsl_matrix* A, gsl_vector* b, unsigned int* basis);
void problem_free(problem_t* pp);

solution_t solve(problem_t problem);

#endif