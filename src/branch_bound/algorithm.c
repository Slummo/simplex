#include "branch_bound/algorithm.h"
#include "branch_bound/stack.h"
#include "simplex/primal.h"
#include "simplex/dual.h"

#include <stdlib.h>

uint32_t solve_relaxation(solve_fn solver, uint32_t is_max, bb_node_t* node_ptr, int32_t* N, solution_t* solution_ptr,
                          uint32_t* iter_n_ptr) {
    return (solver)(node_ptr->state.n, node_ptr->state.m, is_max, &node_ptr->c_view.vector, &node_ptr->A_view.matrix,
                    &node_ptr->b_view.vector, node_ptr->B_view, N, solution_ptr, iter_n_ptr);
}

// Choses a non-integer variable to start branching from.
// Returns -2 on error, -1 if the solution contains only
// integers or the index of the first non-integer
// variable on success
int32_t select_branch_var(const var_arr_t* var_arr_ptr, const solution_t* current_sol_ptr) {
    if (!var_arr_ptr || !current_sol_ptr) {
        return -2;
    }

    for (uint32_t i = 0; i < var_arr_length(var_arr_ptr); i++) {
        const variable_t* v = var_arr_get(var_arr_ptr, i);
        if (variable_is_integer(v) && !solution_var_is_integer(current_sol_ptr, i)) {
            return (int32_t)i;
        }
    }

    return -1;
}

// Branch and bound method on linear problem p
uint32_t branch_and_bound(problem_t* problem_ptr, solution_t* solution_ptr) {
    if (!problem_ptr || !solution_ptr) {
        return 0;
    }

    pstack_t stack = {0};
    if (!pstack_init(&stack)) {
        return 0;
    }

    bb_arena_t arena = {0};
    if (!bb_arena_init(&arena, MAX_N, MAX_M)) {
        pstack_free(&stack);
        return 0;
    }

    bb_arena_copy_problem(&arena, problem_ptr);

    solution_t* best = NULL;
    var_arr_t var_arr = {0};
    if (!var_arr_duplicate(problem_var_arr(problem_ptr), &var_arr)) {
        goto fail;
    }

    bb_node_t root = {0};
    bb_node_init_root(&root, problem_n(problem_ptr), problem_m(problem_ptr), &arena);

    if (!pstack_push(&stack, &root)) {
        goto fail;
    }

    uint32_t is_max = problem_is_max(problem_ptr);
    int32_t* N = problem_N_mut(problem_ptr);

    uint32_t is_root = 1;
    while (!pstack_empty(&stack)) {
        bb_node_t* current_node = pstack_pop(&stack);
        if (!current_node) {
            continue;
        }

        solution_t current_solution = {0};
        uint32_t iter_n = 0;
        solve_fn solver = is_root ? simplex_primal : simplex_dual;

        if (!solve_relaxation(solver, is_max, current_node, N, &current_solution, &iter_n)) {
            goto fail;
        }

        if (is_root) {
            is_root = 0;
        }

        if (solution_is_unbounded(&current_solution)) {
            solution_free(&current_solution);
            continue;
        }

        int32_t branch_var = select_branch_var(&var_arr, &current_solution);
        if (branch_var == -2) {
            solution_free(&current_solution);
            continue;
        }

        // Only integer variables
        if (branch_var == -1) {
            if (!best) {
                best = &current_solution;
            } else if (is_max ? solution_z(&current_solution) > solution_z(best)
                              : solution_z(&current_solution) < solution_z(best)) {
                solution_free(best);
                best = &current_solution;
            } else {
                solution_free(&current_solution);
            }

            continue;
        }

        // Branch on non-integer variable
        double bound = gsl_vector_get(solution_x(&current_solution), branch_var);
        if (bb_node_branch(current_node, &arena, branch_var, bound, 'U', &var_arr)) {
            pstack_push(&stack, current_node);
        }

        bb_node_revert_to_parent(current_node, &arena);

        if (bb_node_branch(current_node, &arena, branch_var, bound, 'L', &var_arr)) {
            pstack_push(&stack, current_node);
        }

        solution_free(&current_solution);
    }

    pstack_free(&stack);
    bb_arena_free(&arena);
    var_arr_free(&var_arr);

    if (best) {
        *solution_ptr = *best;
    }

    return best != NULL;

fail:
    pstack_free(&stack);
    bb_arena_free(&arena);
    var_arr_free(&var_arr);
    return 0;
}