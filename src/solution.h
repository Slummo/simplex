#ifndef SOLUTION_H
#define SOLUTION_H

#include <gsl/gsl_vector.h>
#include <stdint.h>

typedef struct solution _solution, *solution_t;
solution_t solution_new(uint32_t n, gsl_vector* x, uint32_t* basis, uint32_t is_unbounded, uint32_t n_iter);
void solution_print(const solution_t s);
void solution_free(solution_t* sp);
gsl_vector* solution_optimal_vec(const solution_t s);
double solution_optimal_value(const solution_t s);
uint32_t* solution_basis(const solution_t s);
uint32_t is_solution_unbounded(const solution_t s);
uint32_t solution_iterations(const solution_t s);
void solution_set_optimal_value(const solution_t s, double optimal_value);
void solution_set_iterations(const solution_t s, uint32_t iterations);

#endif