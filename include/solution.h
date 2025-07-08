#ifndef SOLUTION_H
#define SOLUTION_H

#include <stdint.h>
#include <gsl/gsl_vector.h>

typedef struct solution {
    uint32_t n;             // Number of constraints
    uint32_t m;             // Number of variables
    gsl_vector* x;          // Optimal solution (m)
    double z;               // Optimal value
    uint32_t is_unbounded;  // Boolean value to know if unbounded
    uint32_t pI_iter;       // Number of iterations of PhaseI to find a base
    uint32_t pII_iter;      // Number of iterations of PhaseII to find solution
} solution_t;

uint32_t solution_init(solution_t* solution_ptr, uint32_t n, uint32_t m, uint32_t is_unbounded, uint32_t pI_iter,
                       uint32_t pII_iter);

// Checks if the i-th component of the solution is an integer
uint32_t solution_var_is_integer(const solution_t* solution_ptr, uint32_t i);

// Pretty print
void solution_print(const solution_t* solution_ptr, const char* name);

void solution_free(solution_t* solution_ptr);

/* GETTERS */
const gsl_vector* solution_x(const solution_t* solution_ptr);
gsl_vector* solution_x_mut(solution_t* solution_ptr);
double solution_z(const solution_t* solution_ptr);
uint32_t solution_is_unbounded(const solution_t* solution_ptr);
uint32_t solution_pI_iterations(const solution_t* solution_ptr);
uint32_t solution_pII_iterations(const solution_t* solution_ptr);

/* SETTERS */
uint32_t solution_set_optimal_value(solution_t* solution_ptr, double optimal_value);
uint32_t solution_set_pI_iter(solution_t* solution_ptr, uint32_t pI_iter);

#endif