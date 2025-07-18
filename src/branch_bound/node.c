#include "branch_bound/node.h"

#include <math.h>

void bb_node_state_init(struct bb_node_state* state_ptr, uint32_t n, uint32_t m) {
    state_ptr->n = n;
    state_ptr->m = m;
}

uint32_t bb_node_init_root(bb_node_t* node_ptr, uint32_t n, uint32_t m, const bb_arena_t* arena_ptr) {
    if (!node_ptr || !arena_ptr) {
        return 0;
    }

    bb_node_state_init(&node_ptr->state, n, m);
    node_ptr->c_view = bb_arena_get_c_view(arena_ptr, n, m);
    node_ptr->A_view = bb_arena_get_A_view(arena_ptr, n, m);
    node_ptr->b_view = bb_arena_get_b_view(arena_ptr, n);
    node_ptr->B_view = bb_arena_get_B_view(arena_ptr, n);
    node_ptr->parent_state = node_ptr->state;

    return 1;
}

uint32_t bb_node_init(bb_node_t* node_ptr, uint32_t n, uint32_t m, const bb_arena_t* arena_ptr,
                      struct bb_node_state parent_state) {
    if (!node_ptr || !arena_ptr) {
        return 0;
    }

    bb_node_state_init(&node_ptr->state, n, m);
    node_ptr->c_view = bb_arena_get_c_view(arena_ptr, n, m);
    node_ptr->A_view = bb_arena_get_A_view(arena_ptr, n, m);
    node_ptr->b_view = bb_arena_get_b_view(arena_ptr, n);
    node_ptr->B_view = bb_arena_get_B_view(arena_ptr, n);
    node_ptr->parent_state = parent_state;

    return 1;
}

// Adds a bound on the variable with index branch_var_index.
// Direction can either be 'L' (Lower bound) or 'U' (Upper bound)
uint32_t bb_node_branch(bb_node_t* node_ptr, const bb_arena_t* arena_ptr, int32_t branch_var_index, double bound,
                        char direction, var_arr_t* var_arr_ptr) {
    if ((direction != 'U' && direction != 'L') || !node_ptr || !arena_ptr || !var_arr_ptr || branch_var_index < 0) {
        fprintf(stderr, "NULL parameters or unknown direction in bb_node_branch\n");
        return 0;
    }

    struct bb_node_state parent_state = {.n = node_ptr->state.n, .m = node_ptr->state.m};
    uint32_t n = parent_state.n;
    uint32_t m = parent_state.m;

    if (!bb_node_init(node_ptr, n + 1, m + 1, arena_ptr, parent_state)) {
        return 0;
    }

    // Direction == 'U' => x[branch_var_index] <= floor(bound)
    // Direction == 'L' => x[branch_var_index] >= ceil(bound)
    bound = direction == 'U' ? floor(bound) : ceil(bound);

    gsl_matrix* A = &node_ptr->A_view.matrix;
    gsl_vector* b = &node_ptr->b_view.vector;
    int32_t* B = node_ptr->B_view;

    // Add bound value
    gsl_vector_set(b, n, bound);

    // Set the branch variable and the new slack/surplus variable
    gsl_matrix_set(A, n, (uint32_t)branch_var_index, 1.0);
    gsl_matrix_set(A, n, m, direction == 'U' ? 1.0 : -1.0);

    // Add the new slack/surplus variable to the base
    B[n] = (int32_t)m;

    // Push the new variable to the array
    variable_t v;
    if (!variable_init_real_positive(&v, 10e9) || !var_arr_push(var_arr_ptr, &v)) {
        return 0;
    }

    return 1;
}

uint32_t bb_node_revert_to_parent(bb_node_t* node_ptr, const bb_arena_t* arena_ptr) {
    struct bb_node_state parent_state = node_ptr->parent_state;
    return bb_node_init(node_ptr, parent_state.n, parent_state.m, arena_ptr, parent_state);
}