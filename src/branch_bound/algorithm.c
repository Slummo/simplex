#include "branch_bound/algorithm.h"
#include "simplex/primal.h"
#include "simplex/dual.h"

#include <stdlib.h>
#include <string.h>

uint32_t init(const problem_t* problem_ptr, pstack_t* stack_ptr, bb_arena_t* arena_ptr, var_arr_t* var_arr_ptr) {
    if (!pstack_init(stack_ptr)) {
        return 0;
    }

    if (!bb_arena_init(arena_ptr, MAX_N, MAX_M)) {
        pstack_free(stack_ptr);
        return 0;
    }

    bb_arena_copy_problem(arena_ptr, problem_ptr);

    if (!var_arr_init(var_arr_ptr, MAX_M)) {
        pstack_free(stack_ptr);
        bb_arena_free(arena_ptr);
        return 0;
    }

    const var_arr_t* var_arr_og = problem_var_arr(problem_ptr);

    memcpy(var_arr_ptr->data, var_arr_og->data, sizeof(variable_t) * var_arr_og->length);
    var_arr_ptr->length = var_arr_og->length;

    return 1;
}

uint32_t solve_relaxation(solve_fn solver, uint32_t is_max, bb_node_t* node_ptr, int32_t* N, solution_t* solution_ptr,
                          uint32_t* iter_n_ptr) {
    return (solver)(node_ptr->state.n, node_ptr->state.m, is_max, &node_ptr->c_view.vector, &node_ptr->A_view.matrix,
                    &node_ptr->b_view.vector, node_ptr->B_view, N, solution_ptr, iter_n_ptr);
}

// Choses a non-integer variable to start branching from.
// Returns -1 if the solution contains only
// integers or the index of the first non-integer
// variable on success
int32_t select_branch_var(const var_arr_t* var_arr_ptr, const solution_t* current_sol_ptr) {
    for (uint32_t i = 0; i < var_arr_length(var_arr_ptr); i++) {
        const variable_t* v = var_arr_get(var_arr_ptr, i);
        if (variable_is_integer(v) && !solution_var_is_integer(current_sol_ptr, i)) {
            return (int32_t)i;
        }
    }

    return -1;
}

uint32_t update(solution_t* best_solution, solution_t* current_solution, bb_node_t* current_node, pstack_t* stack_ptr) {
    double current_z = solution_z(current_solution);
    if (solution_is_unbounded(current_solution) || current_z <= solution_z(best_solution)) {
        // Prune the node
        solution_free(current_solution);
        return 1;
    }

    if (!pstack_push(stack_ptr, *current_node)) {
        return 0;
    }

    if (solution_is_integer(current_solution)) {
        solution_set_x(best_solution, (gsl_vector*)solution_x(current_solution));
        solution_set_z(best_solution, current_z);
    }

    return 1;
}

// Branch and bound method on linear problem p
uint32_t branch_and_bound(problem_t* problem_ptr, solution_t* solution_ptr) {
    if (!problem_ptr || !solution_ptr) {
        return 0;
    }

    pstack_t stack = {0};
    bb_arena_t arena = {0};
    var_arr_t var_arr = {0};
    if (!init(problem_ptr, &stack, &arena, &var_arr)) {
        return 0;
    }

    uint32_t ret = 1;
    uint32_t is_max = problem_is_max(problem_ptr);
    int32_t* N = problem_N_mut(problem_ptr);

    // Start of the algorithm
    solution_t best = {0};

    // Solve root relaxation
    bb_node_t root = {0};
    bb_node_init_root(&root, problem_n(problem_ptr), problem_m(problem_ptr), &arena);
    uint32_t iter_n = 0;
    if (!solve_relaxation(simplex_primal, is_max, &root, N, &best, &iter_n)) {
        goto fail;
    }

    // If the solution of the root relaxation is already integer
    // return it
    solution_print(&best, "Root solution");
    if (solution_is_integer(&best)) {
        goto cleanup;
    }

    // Else push the root to the stack
    if (!pstack_push(&stack, root)) {
        goto fail;
    }
    best.z = -1e20;

    while (!pstack_empty(&stack)) {
        bb_node_t current_node = {0};
        pstack_pop(&stack, &current_node);

        // Branch, solve relaxations and push nodes into the stack
        int32_t branch_var = select_branch_var(&var_arr, &best);
        if (branch_var == -1) {
            goto fail;
        }
        double bound = gsl_vector_get(solution_x(&best), branch_var);

        // Left branch
        if (!bb_node_branch(&current_node, &arena, branch_var, bound, 'U', &var_arr)) {
            goto fail;
        }
        solution_t left_solution = {0};
        uint32_t left_iter_n = 0;
        if (!solve_relaxation(simplex_dual, is_max, &current_node, N, &left_solution, &left_iter_n)) {
            goto fail;
        }
        if (!update(&best, &left_solution, &current_node, &stack)) {
            goto fail;
        }
        bb_node_revert_to_parent(&current_node, &arena);

        // Right branch
        if (!bb_node_branch(&current_node, &arena, branch_var, bound, 'U', &var_arr)) {
            goto fail;
        }
        solution_t right_solution = {0};
        uint32_t right_iter_n = 0;
        if (!solve_relaxation(simplex_dual, is_max, &current_node, N, &right_solution, &right_iter_n)) {
            goto fail;
        }
        if (!update(&best, &right_solution, &current_node, &stack)) {
            goto fail;
        }
        bb_node_revert_to_parent(&current_node, &arena);
    }

    goto cleanup;

fail:
    ret = 0;

cleanup:
    pstack_free(&stack);
    bb_arena_free(&arena);
    var_arr_free(&var_arr);
    *solution_ptr = best;
    return ret;
}