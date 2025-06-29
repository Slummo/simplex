#ifndef PROBLEM_H
#define PROBLEM_H

#define MAX_ROWS 100

#include "solution.h"
#include "variable.h"
#include <stdint.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

typedef struct problem problem_t;

problem_t* problem_new(uint32_t n, uint32_t m, uint32_t is_max, const gsl_vector* c_raw, const gsl_matrix* A_raw,
                       const gsl_vector* b_raw, const int32_t* basis_raw, uint32_t pI_iter, variable_t** variables_raw);

problem_t* problem_new2(uint32_t n, uint32_t m, uint32_t is_max, const gsl_vector* c_raw, const gsl_matrix* A_raw,
                        const gsl_vector* b_raw, variable_t** variables_raw);

problem_t* problem_duplicate(const problem_t* p);

uint32_t problem_is_milp(const problem_t* p);

// Pretty print
void problem_print(const problem_t* p, const char* name);

void problem_free(problem_t** pp);

/* GETTERS */
uint32_t problem_n(const problem_t* p);
uint32_t problem_m(const problem_t* p);
uint32_t problem_is_max(const problem_t* p);
const gsl_vector* problem_c(const problem_t* p);
const gsl_matrix* problem_A(const problem_t* p);
const gsl_vector* problem_b(const problem_t* p);
const int32_t* problem_basis(const problem_t* p);
int32_t* problem_basis_mut(const problem_t* p);
const int32_t* problem_nonbasis(const problem_t* p);
int32_t* problem_nonbasis_mut(const problem_t* p);
uint32_t problem_pI_iter(const problem_t* p);
const varr_t* problem_varr(const problem_t* p);

#endif