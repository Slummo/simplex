#ifndef PROBLEM_H
#define PROBLEM_H

#define MAX_ROWS 100

#include "vector.h"
#include "matrix.h"
#include "variable.h"
#include "solution.h"
#include <stdint.h>

typedef struct problem {
    uint32_t n;         // Number of constraints
    uint32_t m;         // Number of variables
    uint32_t is_max;    // Boolean value to know if its a maximization problem
    vector_t c;         // Reduced costs (m)
    matrix_t A;         // Constraints matrix (n x m)
    vector_t b;         // RHS (n)
    int32_t* B;         // Indices of basic variables (size n)
    int32_t* N;         // Indices of nonbasic variables (size m-n)
    uint32_t pI_iter;   // Number of iterations to find base with PhaseI
    var_arr_t var_arr;  // Array of variables
} problem_t;

// uint32_t problem_init(problem_t* problem_ptr, uint32_t n, uint32_t m, uint32_t is_max, gsl_vector* c, gsl_matrix* A,
//                       gsl_vector* b, int32_t* B, uint32_t pI_iter, variable_t* variables);

// uint32_t problem_init2(problem_t* problem_ptr, uint32_t n, uint32_t m, uint32_t is_max, gsl_vector* c, gsl_matrix* A,
//                        gsl_vector* b, variable_t* variables);

uint32_t problem_from_model(problem_t* problem_ptr, FILE* stream);

uint32_t problem_is_milp(const problem_t* problem_ptr);

// Pretty print
void problem_print(const problem_t* problem_ptr, const char* name);

/// @brief Checks if a problem has a feasible base
/// @param problem_ptr A const pointer to the problem to checl
/// @param B A dinamically allocated basis indices array
/// @return 1 if the problem has a feasible base else 0
uint32_t problem_has_feasible_base(const problem_t* problem_ptr, int32_t* B);

uint32_t solve_with_simplex(problem_t* problem_ptr, solution_t* solution_ptr);

void problem_free(problem_t* problem_ptr);

/* GETTERS */
uint32_t problem_n(const problem_t* problem_ptr);
uint32_t problem_m(const problem_t* problem_ptr);
uint32_t problem_is_max(const problem_t* problem_ptr);
const vector_t* problem_c(const problem_t* problem_ptr);
vector_t* problem_c_mut(problem_t* problem_ptr);
uint32_t problem_c_as_gsl_view(problem_t* problem_ptr, uint32_t offset, uint32_t length, gsl_vector_view* view_ptr);
const matrix_t* problem_A(const problem_t* problem_ptr);
matrix_t* problem_A_mut(problem_t* problem_ptr);
uint32_t problem_A_as_gsl_view(problem_t* problem_ptr, uint32_t row_offset, uint32_t col_offset, uint32_t rows,
                               uint32_t cols, gsl_matrix_view* view_ptr);
const vector_t* problem_b(const problem_t* problem_ptr);
vector_t* problem_b_mut(problem_t* problem_ptr);
uint32_t problem_b_as_gsl_view(problem_t* problem_ptr, uint32_t offset, uint32_t length, gsl_vector_view* view_ptr);
const int32_t* problem_B(const problem_t* problem_ptr);
int32_t* problem_B_mut(problem_t* problem_ptr);
const int32_t* problem_N(const problem_t* problem_ptr);
int32_t* problem_N_mut(problem_t* problem_ptr);
uint32_t problem_pI_iter(const problem_t* problem_ptr);
const var_arr_t* problem_var_arr(const problem_t* problem_ptr);
var_arr_t* problem_var_arr_mut(problem_t* problem_ptr);

/* SETTERS */
uint32_t problem_set_pI_iter(problem_t* problem_ptr, uint32_t pI_iter);

#endif