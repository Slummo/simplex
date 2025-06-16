#ifndef PROBLEM_H
#define PROBLEM_H

#define MAX_LINE 1024
#define MAX_ROWS 100

#include <stdio.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>

typedef struct problem _problem, *problem_t;
problem_t problem_new(int n, int m, gsl_vector* c, gsl_matrix* A, gsl_vector* b, int* basis);
void problem_free(problem_t* pp);

typedef struct solution _solution, *solution_t;
void solution_print(solution_t s);
void solution_free(solution_t* sp);

solution_t solve(problem_t problem);

#endif