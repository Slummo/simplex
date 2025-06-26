#ifndef PROBLEM_H
#define PROBLEM_H

#define MAX_ROWS 100

#include "solution.h"
#include <stdint.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

typedef struct problem _problem, *problem_t;

// Duplicates each mallocable param
problem_t problem_new(uint32_t n, uint32_t m, uint32_t is_max, const gsl_vector* c, const gsl_matrix* A,
                      const gsl_vector* b, const int32_t* basis, uint32_t pI_iter, const uint32_t* is_integer);

problem_t problem_new2(uint32_t n, uint32_t m, uint32_t is_max, const gsl_vector* c, const gsl_matrix* A,
                       const gsl_vector* b, const uint32_t* is_integer);

problem_t problem_duplicate(const problem_t p);

uint32_t problem_is_milp(const problem_t p);

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
uint32_t problem_pI_iter(const problem_t p);
const uint32_t* problem_integers(const problem_t p);

/* SETTERS */
void problem_set_basis_nonbasis(const problem_t p, int32_t* basis);

#endif