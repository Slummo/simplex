#include "problem.h"
#include "utils.h"
#include "simplex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

uint32_t problem_from_model(problem_t* problem_ptr, FILE* stream) {
    if (!problem_ptr || !stream) {
        return 0;
    }

    uint32_t n;
    uint32_t m;
    uint32_t is_max;
    vector_t c;
    matrix_t A;
    vector_t b;
    var_arr_t var_arr;
    int32_t* B = NULL;
    int32_t* N = NULL;

    if (stream == stdin) {
        printf("n: ");
        fflush(stdout);
    }

    if (fscanf(stream, "%u", &n) != 1 || n > MAX_ROWS) {
        fprintf(stderr, "Invalid or too many rows (max %u)\n", MAX_ROWS);
        goto fail;
    }

    if (stream == stdin) {
        printf("m: ");
        fflush(stdout);
    }
    if (fscanf(stream, "%u", &m) != 1) {
        fprintf(stderr, "Failed to read m\n");
        goto fail;
    }
    fgetc(stream);  // consume newline

    if (stream == stdin) {
        printf("is_max: ");
        fflush(stdout);
    }
    if (fscanf(stream, "%u", &is_max) != 1) {
        fprintf(stderr, "Failed to read is_max\n");
        goto fail;
    }

    // Augmented capacity for phaseI
    uint32_t variables_num = m + n;
    uint32_t constraints_num = n;

    if (!vector_from_stream(&c, stream, "c", variables_num, m)) {
        goto fail;
    }

    if (!matrix_from_stream(&A, stream, "A", constraints_num, variables_num, n, m)) {
        goto fail;
    }

    if (!vector_from_stream(&b, stream, "b", constraints_num, n)) {
        goto fail;
    }

    if (!var_arr_from_stream(&var_arr, stream, variables_num, m)) {
        goto fail;
    }

    problem_ptr->n = n;
    problem_ptr->m = m;
    problem_ptr->is_max = is_max;
    problem_ptr->c = c;
    problem_ptr->A = A;
    problem_ptr->b = b;
    problem_ptr->var_arr = var_arr;

    // Find a base with phaseI
    B = simplex_phaseI(problem_ptr);
    if (!B) {
        goto fail;
    }

    N = calculate_nonbasis(B, n, m);
    if (!N) {
        goto fail;
    }

    problem_ptr->B = B;
    problem_ptr->N = N;

    if (stream != stdin) {
        fclose(stream);
    }

    return 1;

fail:
    vector_free(&c);
    matrix_free(&A);
    vector_free(&b);
    var_arr_free(&var_arr);
    free(B);
    free(N);
    if (stream != stdin) {
        fclose(stream);
    }
    return 0;
}

uint32_t problem_is_milp(const problem_t* problem_ptr) {
    if (!problem_ptr) {
        return 0;
    }

    const var_arr_t* var_arr = &problem_ptr->var_arr;
    for (uint32_t i = 0; i < problem_ptr->m; i++) {
        if (var_arr_get(var_arr, i)->type == VAR_INTEGER) {
            return 1;
        }
    }

    return 0;
}

/// @brief Checks if a problem has a feasible base. If not the B array is memset to 0
/// @param problem_ptr A const pointer to the problem to checl
/// @param B A dinamically allocated basis indices array
/// @return 1 if the problem has a feasible base else 0
uint32_t problem_has_feasible_base(const problem_t* problem_ptr, int32_t* B) {
    if (!problem_ptr || !B) {
        return 0;
    }

    uint32_t n = problem_n(problem_ptr);
    uint32_t m = problem_m(problem_ptr);

    // Initialize basis indices to -1
    for (uint32_t i = 0; i < n; i++) {
        B[i] = -1;
    }

    // Check for identity
    const matrix_t* A = problem_A(problem_ptr);
    for (uint32_t j = 0; j < m; j++) {
        int32_t pivot_row = -1;
        uint32_t valid = 1;
        for (uint32_t i = 0; i < n; i++) {
            double aij = matrix_get(A, i, j);
            if (aij == 1.0) {
                if (pivot_row == -1) {
                    pivot_row = (int32_t)i;
                } else {
                    valid = 0;
                    break;
                }
            } else if (aij != 0.0) {
                valid = 0;
                break;
            }
        }
        if (valid && pivot_row != -1 && B[pivot_row] == -1) {
            B[pivot_row] = j;
        }
    }

    // Check if base is feasible
    const vector_t* b = problem_b(problem_ptr);
    uint32_t is_feasible = 1;
    for (uint32_t i = 0; is_feasible && i < n; i++) {
        if (B[i] == -1 || vector_get(b, i) < 0) {
            is_feasible = 0;
        }
    }

    if (!is_feasible) {
        memset(B, 0, sizeof(int32_t) * n);
    }

    return is_feasible;
}

uint32_t solve_with_simplex(problem_t* problem_ptr, solution_t* solution_ptr) {
    uint32_t n = problem_n(problem_ptr);
    uint32_t m = problem_m(problem_ptr);
    uint32_t is_max = problem_is_max(problem_ptr);

    // Get views
    gsl_vector_view c_view;
    if (!problem_c_as_gsl_view(problem_ptr, 0, m, &c_view)) {
        return 0;
    }
    gsl_vector* c_gsl = &c_view.vector;

    gsl_matrix_view A_view;
    if (!problem_A_as_gsl_view(problem_ptr, 0, 0, n, m, &A_view)) {
        return 0;
    }
    gsl_matrix* A_gsl = &A_view.matrix;

    gsl_vector_view b_view;
    if (!problem_b_as_gsl_view(problem_ptr, 0, n, &b_view)) {
        return 0;
    }
    gsl_vector* b_gsl = &b_view.vector;

    return simplex_phaseII(n, m, is_max, c_gsl, A_gsl, b_gsl, problem_B_mut(problem_ptr), problem_N_mut(problem_ptr),
                           solution_ptr);
}

#define TERM_WIDTH 8

void print_coefficient(double coef, uint32_t var_idx, int is_first) {
    char buf[TERM_WIDTH + 1];

    if (fabs(coef) < 1e-9) {
        // Print spaces for alignment
        snprintf(buf, sizeof(buf), "%*s", TERM_WIDTH, "");
    } else {
        char sign = coef < 0 ? '-' : (is_first ? ' ' : '+');
        double abs_val = fabs(coef);

        // Check if it's 1 or -1
        if (fabs(abs_val - 1.0) < 1e-9) {
            snprintf(buf, sizeof(buf), "%c x%-2u", sign, var_idx + 1);  // omit the 1
        } else if (fabs(abs_val - (int32_t)abs_val) < 1e-9) {
            snprintf(buf, sizeof(buf), "%c%2dx%-2u", sign, (int32_t)abs_val, var_idx + 1);
        } else {
            snprintf(buf, sizeof(buf), "%c%.1lfx%-2u", sign, abs_val, var_idx + 1);
        }
    }

    printf("%*s", TERM_WIDTH, buf);
}

// Pretty print
void problem_print(const problem_t* problem_ptr, const char* name) {
    if (!problem_ptr || !name) {
        return;
    }

    printf("\n================== %s ==================\n", name);
    printf("%s z = ", problem_ptr->is_max ? "max" : "min");

    // Objective function
    const vector_t* c = &problem_ptr->c;
    for (uint32_t i = 0; i < problem_ptr->m; i++) {
        double ci = vector_get(c, i);
        print_coefficient(ci, i, i == 0);
    }
    printf("\n\nconstraints:\n");

    // Constraints
    const matrix_t* A = &problem_ptr->A;
    const vector_t* b = &problem_ptr->b;
    for (uint32_t i = 0; i < problem_ptr->n; i++) {
        printf("\t");
        for (uint32_t j = 0; j < problem_ptr->m; j++) {
            double aij = matrix_get(A, i, j);
            print_coefficient(aij, j, j == 0);
        }

        double bi = vector_get(b, i);
        printf(" = ");
        if (fabs(bi - (int32_t)bi) < 1e-9) {
            printf("%d", (int32_t)bi);
        } else {
            printf("%.3lf", bi);
        }
        printf("\n");
    }

    printf("\nvariables:\n");
    const var_arr_t* var_arr = &problem_ptr->var_arr;
    for (uint32_t i = 0; i < problem_ptr->m; i++) {
        printf("\t");
        variable_print(var_arr_get(var_arr, i));
    }
    puts("");
}

void problem_free(problem_t* problem_ptr) {
    if (!problem_ptr) {
        return;
    }

    vector_free(&problem_ptr->c);
    matrix_free(&problem_ptr->A);
    vector_free(&problem_ptr->b);
    free(problem_ptr->B);
    free(problem_ptr->N);
    var_arr_free(&problem_ptr->var_arr);
}

/* GETTERS */

uint32_t problem_n(const problem_t* problem_ptr) {
    return problem_ptr ? problem_ptr->n : 0;
}

uint32_t problem_m(const problem_t* problem_ptr) {
    return problem_ptr ? problem_ptr->m : 0;
}

uint32_t problem_is_max(const problem_t* problem_ptr) {
    return problem_ptr ? problem_ptr->is_max : 0;
}

const vector_t* problem_c(const problem_t* problem_ptr) {
    return problem_ptr ? &problem_ptr->c : NULL;
}

vector_t* problem_c_mut(problem_t* problem_ptr) {
    return problem_ptr ? &problem_ptr->c : NULL;
}

uint32_t problem_c_as_gsl_view(problem_t* problem_ptr, uint32_t offset, uint32_t length, gsl_vector_view* view_ptr) {
    return vector_as_gsl_view(&problem_ptr->c, offset, length, view_ptr);
}

const matrix_t* problem_A(const problem_t* problem_ptr) {
    return problem_ptr ? &problem_ptr->A : NULL;
}

matrix_t* problem_A_mut(problem_t* problem_ptr) {
    return problem_ptr ? &problem_ptr->A : NULL;
}

uint32_t problem_A_as_gsl_view(problem_t* problem_ptr, uint32_t row_offset, uint32_t col_offset, uint32_t rows,
                               uint32_t cols, gsl_matrix_view* view_ptr) {
    return matrix_as_gsl_view(&problem_ptr->A, row_offset, col_offset, rows, cols, view_ptr);
}

const vector_t* problem_b(const problem_t* problem_ptr) {
    return problem_ptr ? &problem_ptr->b : NULL;
}

vector_t* problem_b_mut(problem_t* problem_ptr) {
    return problem_ptr ? &problem_ptr->b : NULL;
}

uint32_t problem_b_as_gsl_view(problem_t* problem_ptr, uint32_t offset, uint32_t length, gsl_vector_view* view_ptr) {
    return vector_as_gsl_view(&problem_ptr->b, offset, length, view_ptr);
}

const int32_t* problem_B(const problem_t* problem_ptr) {
    return problem_ptr ? problem_ptr->B : NULL;
}

int32_t* problem_B_mut(problem_t* problem_ptr) {
    return problem_ptr ? problem_ptr->B : NULL;
}

const int32_t* problem_N(const problem_t* problem_ptr) {
    return problem_ptr ? problem_ptr->N : NULL;
}

int32_t* problem_N_mut(problem_t* problem_ptr) {
    return problem_ptr ? problem_ptr->N : NULL;
}

uint32_t problem_pI_iter(const problem_t* problem_ptr) {
    return problem_ptr ? problem_ptr->pI_iter : 0;
}

const var_arr_t* problem_var_arr(const problem_t* problem_ptr) {
    return problem_ptr ? &problem_ptr->var_arr : NULL;
}

var_arr_t* problem_var_arr_mut(problem_t* problem_ptr) {
    return problem_ptr ? &problem_ptr->var_arr : NULL;
}

/* SETTERS */
uint32_t problem_set_pI_iter(problem_t* problem_ptr, uint32_t pI_iter) {
    if (!problem_ptr) {
        return 0;
    }

    problem_ptr->pI_iter = pI_iter;

    return 1;
}