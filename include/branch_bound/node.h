#ifndef BB_NODE_H
#define BB_NODE_H

#include "branch_bound/arena.h"
#include "variable.h"
#include "problem.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

struct bb_node_state {
    uint32_t n;
    uint32_t m;
};

void bb_node_state_init(struct bb_node_state* state_ptr, uint32_t n, uint32_t m);

typedef struct bb_node {
    struct bb_node_state state;
    gsl_vector_view c_view;
    gsl_matrix_view A_view;
    gsl_vector_view b_view;
    int32_t* B_view;
    struct bb_node_state parent_state;
} bb_node_t;

uint32_t bb_node_init_root(bb_node_t* node_ptr, uint32_t n, uint32_t m, const bb_arena_t* arena_ptr);

uint32_t bb_node_init(bb_node_t* node_ptr, uint32_t n, uint32_t m, const bb_arena_t* arena_ptr,
                      struct bb_node_state parent_state);

// Adds a bound on the variable with index branch_var_index.
// Direction can either be 'L' (Lower bound) or 'U' (Upper bound)
uint32_t bb_node_branch(bb_node_t* node_ptr, const bb_arena_t* arena_ptr, int32_t branch_var_index, double bound,
                        char direction, var_arr_t* var_arr_ptr);

uint32_t bb_node_revert_to_parent(bb_node_t* node_ptr, const bb_arena_t* arena_ptr);

#endif