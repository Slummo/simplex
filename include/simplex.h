#ifndef SIMPLEX_H
#define SIMPLEX_H

#include "problem.h"
#include "solution.h"
#include "variable.h"
#include <stdint.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

// Find problem basis indices with Phase 1 method
int32_t* simplex_phaseI(problem_t* problem_ptr, uint32_t* iter_n_ptr);

void extract_basic_objects(uint32_t n, uint32_t is_max, int32_t* B, const gsl_vector* c, const gsl_matrix* A,
                           gsl_vector* cB, gsl_matrix* AB);
gsl_matrix* inverse(const gsl_matrix* base, size_t size);
void compute_basic_solution(const gsl_matrix* AB_inv, const gsl_vector* b, gsl_vector* xB);
void compute_reduced_costs(uint32_t n, uint32_t m, uint32_t is_max, int32_t* N, const gsl_vector* c, gsl_vector* cB,
                           gsl_vector* cN, const gsl_matrix* A, const gsl_matrix* AB_inv, gsl_vector* r);
uint32_t extract_column(const gsl_matrix* m, uint32_t j, gsl_vector* col);
uint32_t extract_row(const gsl_matrix* m, uint32_t i, gsl_vector* row);
void pivot(int32_t entering, int32_t leaving, int32_t* B, int32_t* N);
void extract_optimal(uint32_t n, uint32_t is_max, int32_t* B, gsl_vector* xB, const gsl_vector* c,
                     solution_t* solution_ptr);

uint32_t simplex_primal(uint32_t n, uint32_t m, uint32_t is_max, const gsl_vector* c, const gsl_matrix* A,
                        const gsl_vector* b, int32_t* B, int32_t* N, solution_t* solution_ptr, uint32_t* iter_n_ptr);

uint32_t simplex_dual(uint32_t n, uint32_t m, uint32_t is_max, const gsl_vector* c, const gsl_matrix* A,
                      const gsl_vector* b, int32_t* B, int32_t* N, solution_t* solution_ptr, uint32_t* iter_n_ptr);

#endif