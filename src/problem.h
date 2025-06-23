#ifndef PROBLEM_H
#define PROBLEM_H

#define MAX_ROWS 100

#include "solution.h"
#include <stdint.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

typedef struct problem _problem, *problem_t;

// Frees arguments on error
problem_t problem_new(uint32_t n, uint32_t m, uint32_t is_max, gsl_vector* c, gsl_matrix* A, gsl_vector* b,
                      int32_t* basis);

problem_t problem_duplicate(const problem_t p);

// Pretty print
void problem_print(const problem_t p, const char* name);

void problem_free(problem_t* pp);

/* GETTERS */
uint32_t problem_n(const problem_t p);
uint32_t problem_m(const problem_t p);
uint32_t problem_is_max(const problem_t p);
const gsl_vector* problem_c(const problem_t p);
const gsl_matrix* problem_A(const problem_t p);
const gsl_vector* problem_b(const problem_t p);
const int32_t* problem_basis(const problem_t p);
int32_t* problem_basis_mut(const problem_t p);
const int32_t* problem_nonbasis(const problem_t p);
int32_t* problem_nonbasis_mut(const problem_t p);

#endif