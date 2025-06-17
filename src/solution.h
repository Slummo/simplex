#ifndef SOLUTION_H
#define SOLUTION_H

#include <gsl/gsl_vector.h>

typedef struct solution _solution, *solution_t;
solution_t solution_new(int n, gsl_vector* x, double z, int is_unbounded, int n_iter);
void solution_print(solution_t s);
void solution_free(solution_t* sp);
gsl_vector* solution_optimal_vec_ptr(const solution_t s);
double* solution_optimal_value_ptr(const solution_t s);
int is_solution_unbounded(const solution_t s);

#endif