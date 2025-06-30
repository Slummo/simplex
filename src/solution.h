#ifndef SOLUTION_H
#define SOLUTION_H

#include <stdint.h>
#include <gsl/gsl_vector.h>

typedef struct solution solution_t;

// Duplicates each mallocable param
solution_t* solution_new(uint32_t n, uint32_t m, uint32_t is_unbounded, uint32_t pI_iter, uint32_t pII_iter);

// Checks if the i-th component of the solution is an integer
uint32_t solution_var_is_integer(const solution_t* s, uint32_t i);

// Pretty print
void solution_print(const solution_t* s, const char* name);

void solution_free(solution_t** sp);

/* GETTERS */
const gsl_vector* solution_x(const solution_t* s);
gsl_vector* solution_x_mut(const solution_t* s);
double solution_z(const solution_t* s);
uint32_t solution_is_unbounded(const solution_t* s);
uint32_t solution_pI_iterations(const solution_t* s);
uint32_t solution_pII_iterations(const solution_t* s);

/* SETTERS */
void solution_set_optimal_value(solution_t* s, double optimal_value);

#endif