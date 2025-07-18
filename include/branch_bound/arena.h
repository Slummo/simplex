#ifndef BB_ARENA_H
#define BB_ARENA_H

#include "problem.h"

#include <stdio.h>
#include <stdint.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef struct bb_arena {
    double* data;
    size_t max_n;
    size_t max_m;
    double* c_base;
    double* A_base;
    double* b_base;

    int32_t* B;
} bb_arena_t;

uint32_t bb_arena_init(bb_arena_t* arena_ptr, size_t max_n, size_t max_m);

void bb_arena_copy_problem(bb_arena_t* arena_ptr, const problem_t* problem_ptr);

// Creates view to cost vector c of size (m + n). arena_ptr must not be null
gsl_vector_view bb_arena_get_c_view(const bb_arena_t* arena_ptr, size_t m, size_t n);

// Creates view to constraint matrix A of size (n x (m + n)). arena_ptr must not be null
gsl_matrix_view bb_arena_get_A_view(const bb_arena_t* arena_ptr, size_t n, size_t m);

// Creates view to RHS vector b of size n. arena_ptr must not be null
gsl_vector_view bb_arena_get_b_view(const bb_arena_t* arena_ptr, size_t n);

int32_t* bb_arena_get_B_view(const bb_arena_t* arena_ptr, size_t n);

void bb_arena_free(bb_arena_t* arena_ptr);

#endif