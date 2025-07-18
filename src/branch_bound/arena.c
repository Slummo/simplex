#include "branch_bound/arena.h"

#include <string.h>

uint32_t bb_arena_init(bb_arena_t* arena_ptr, size_t max_n, size_t max_m) {
    if (!arena_ptr) {
        return 0;
    }

    // Account for artificial variables that could be needed
    size_t c_size = max_m + max_n;
    size_t A_size = max_n * (max_m + max_n);
    size_t b_size = max_n;

    size_t total_size = c_size + A_size + b_size;

    arena_ptr->data = (double*)calloc(total_size, sizeof(double));
    if (!arena_ptr->data) {
        return 0;
    }

    arena_ptr->max_n = max_n;
    arena_ptr->max_m = max_m;

    // Memory layout: [c][A][b]
    size_t c_off = 0;
    size_t A_off = c_off + c_size;
    size_t b_off = A_off + A_size;

    arena_ptr->c_base = arena_ptr->data + c_off;
    arena_ptr->A_base = arena_ptr->data + A_off;
    arena_ptr->b_base = arena_ptr->data + b_off;

    arena_ptr->B = (int32_t*)calloc(max_n, sizeof(int32_t));
    if (!arena_ptr->B) {
        free(arena_ptr->data);
        return 0;
    }

    return 1;
}

void bb_arena_copy_problem(bb_arena_t* arena_ptr, const problem_t* problem_ptr) {
    uint32_t n = problem_n(problem_ptr);
    uint32_t m = problem_m(problem_ptr);
    const gsl_vector* c = problem_c(problem_ptr);
    const gsl_matrix* A = problem_A(problem_ptr);
    const gsl_vector* b = problem_b(problem_ptr);
    const int32_t* B = problem_B(problem_ptr);

    memcpy(arena_ptr->c_base, c->data, sizeof(double) * m);
    memcpy(arena_ptr->A_base, A->data, sizeof(double) * n * m);
    memcpy(arena_ptr->b_base, b->data, sizeof(double) * n);

    memcpy(arena_ptr->B, B, sizeof(int32_t) * n);
}

// Creates view to cost vector c of size (m + n)
gsl_vector_view bb_arena_get_c_view(const bb_arena_t* arena_ptr, size_t n, size_t m) {
    return gsl_vector_view_array(arena_ptr->c_base, m + n);
}

// Creates view to constraint matrix A of size (n x (m + n))
gsl_matrix_view bb_arena_get_A_view(const bb_arena_t* arena_ptr, size_t n, size_t m) {
    return gsl_matrix_view_array(arena_ptr->A_base, n, m + n);
}

// Creates view to RHS vector b of size n
gsl_vector_view bb_arena_get_b_view(const bb_arena_t* arena_ptr, size_t n) {
    return gsl_vector_view_array(arena_ptr->b_base, n);
}

int32_t* bb_arena_get_B_view(const bb_arena_t* arena_ptr, size_t n) {
    return arena_ptr->B;
}

void bb_arena_free(bb_arena_t* arena_ptr) {
    if (!arena_ptr) {
        return;
    }

    free(arena_ptr->data);
    arena_ptr->data = NULL;
    arena_ptr->max_n = 0;
    arena_ptr->max_m = 0;
    arena_ptr->c_base = NULL;
    arena_ptr->A_base = NULL;
    arena_ptr->b_base = NULL;
    free(arena_ptr->B);
    arena_ptr->B = NULL;
}