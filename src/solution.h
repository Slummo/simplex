#ifndef SOLUTION_H
#define SOLUTION_H

#include <stdint.h>
#include <gsl/gsl_vector.h>

typedef struct solution _solution, *solution_t;

solution_t solution_new(uint32_t n, gsl_vector* x, int32_t* basis, uint32_t is_unbounded, uint32_t pI_iter,
                        uint32_t pII_iter);

solution_t solution_duplicate(const solution_t s);

// Checks if the i-th component of the solution is an integer
uint32_t solution_var_is_integer(const solution_t s, uint32_t i);

// Pretty print
void solution_print(const solution_t s);

void solution_free(solution_t* sp);

/* GETTERS */
const gsl_vector* solution_x(const solution_t s);
gsl_vector* solution_x_mut(const solution_t s);
double solution_z(const solution_t s);
const int32_t* solution_basis(const solution_t s);
uint32_t is_solution_unbounded(const solution_t s);
uint32_t solution_pI_iterations(const solution_t s);
uint32_t solution_pII_iterations(const solution_t s);

/* SETTERS */
void solution_set_optimal_value(const solution_t s, double optimal_value);

#endif