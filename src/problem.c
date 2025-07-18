#include "problem.h"
#include "utils.h"
#include "simplex/primal.h"
#include "branch_bound/algorithm.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define PRINT_TERM_WIDTH 8

gsl_vector* gsl_vector_from_stream(FILE* stream, char* name, uint32_t capacity, uint32_t size) {
    if (size > capacity) {
        fprintf(stderr, "Requested size exceeds capacity in gsl_vector_from_stream for %s\n", name);
        return NULL;
    }

    gsl_vector* v = gsl_vector_calloc(capacity);
    if (!v) {
        fprintf(stderr, "Failed to allocate gsl_vector %s in gsl_vector_from_stream\n", name);
        return NULL;
    }

    if (stream == stdin) {
        printf("Enter %u values for vector %s: ", size, name);
        fflush(stdout);
    }

    for (uint32_t i = 0; i < size; i++) {
        if (fscanf(stream, "%lf", v->data + i) != 1) {
            fprintf(stderr, "Failed to read element %u of gsl_vector %s\n", i, name);
            gsl_vector_free(v);
            return NULL;
        }
    }

    return v;
}

gsl_matrix* gsl_matrix_from_stream(FILE* stream, char* name, uint32_t row_capacity, uint32_t col_capacity,
                                   uint32_t rows, uint32_t cols) {
    if (rows > row_capacity || cols > col_capacity) {
        fprintf(stderr, "Requested size exceeds capacity in gsl_matrix_from_stream for %s\n", name);
        return NULL;
    }

    gsl_matrix* m = gsl_matrix_calloc(row_capacity, col_capacity);
    if (!m) {
        fprintf(stderr, "Failed to allocate gsl_matrix %s in gsl_matrix_from_stream\n", name);
        return NULL;
    }

    if (stream == stdin) {
        printf("Enter %ux%u = %u values for matrix %s: ", rows, cols, rows * cols, name);
        fflush(stdout);
    }

    for (uint32_t i = 0; i < rows; i++) {
        for (uint32_t j = 0; j < cols; j++) {
            if (fscanf(stream, "%lf", m->data + i * m->tda + j) != 1) {
                fprintf(stderr, "Failed to read element (%u, %u) of gsl_matrix %s\n", i, j, name);
                gsl_matrix_free(m);
                return NULL;
            }
        }
    }

    return m;
}

void problem_make_RHS_positive(uint32_t n, gsl_matrix* A, gsl_vector* b) {
    for (uint32_t i = 0; i < n; i++) {
        double bi = gsl_vector_get(b, i);
        if (bi < 0.0) {
            // Flip bi
            gsl_vector_set(b, i, -bi);

            // Flip row Ai
            gsl_vector_view row = gsl_matrix_row(A, i);
            gsl_vector_scale(&row.vector, -1.0);
        }
    }
}

int32_t* problem_find_primal_base(uint32_t n, uint32_t m, uint32_t is_max, gsl_vector* c, gsl_matrix* A, gsl_vector* b,
                                  var_arr_t* var_arr_ptr, uint32_t* iter_n_ptr) {
    if (!c || !A || !b || !iter_n_ptr) {
        fprintf(stderr, "Some arguments are NULL in problem_find_primal_base\n");
        return NULL;
    }

    int32_t* B = (int32_t*)malloc(sizeof(int32_t) * n);
    if (!B) {
        return NULL;
    }

    // Initialize basis indices to -1
    for (uint32_t i = 0; i < n; i++) {
        B[i] = -1;
    }

    // Check for identity
    uint32_t n_indices_set = 0;
    for (uint32_t j = 0; j < m; j++) {
        int32_t pivot_row = -1;
        uint32_t valid = 1;
        for (uint32_t i = 0; valid && i < n; i++) {
            double aij = gsl_matrix_get(A, i, j);
            if (aij == 1.0) {
                if (pivot_row == -1) {
                    pivot_row = (int32_t)i;
                } else {
                    valid = 0;
                }
            } else if (aij != 0.0) {
                valid = 0;
            }
        }
        if (valid && pivot_row != -1 && B[pivot_row] == -1) {
            B[pivot_row] = j;
            n_indices_set++;
        }
    }

    // If the B array of indices has been filled
    if (n_indices_set == n) {
        return B;
    }

    memset(B, 0, sizeof(int32_t) * n);

    if (!simplex_primal_phaseI(n, m, is_max, c, A, b, B, var_arr_ptr, iter_n_ptr)) {
        fprintf(stderr, "Failed to find initial primal feasible base with PhaseI\n");
    }

    return B;
}

uint32_t problem_from_stream(problem_t* problem_ptr, FILE* stream) {
    if (!problem_ptr || !stream) {
        return 0;
    }

    uint32_t n;
    uint32_t m;
    uint32_t is_max;
    gsl_vector* c = NULL;
    gsl_matrix* A = NULL;
    gsl_vector* b = NULL;
    var_arr_t var_arr = {0};
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

    c = gsl_vector_from_stream(stream, "c", variables_num, m);
    if (!c) {
        goto fail;
    }

    A = gsl_matrix_from_stream(stream, "A", constraints_num, variables_num, n, m);
    if (!A) {
        goto fail;
    }

    b = gsl_vector_from_stream(stream, "b", constraints_num, n);
    if (!b) {
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
    problem_ptr->pI_iter = 0;

    problem_make_RHS_positive(n, A, b);

    B = problem_find_primal_base(n, m, is_max, c, A, b, &var_arr, &problem_ptr->pI_iter);
    if (!B) {
        goto fail;
    }

    N = calculate_nonbasis(B, n, m + n);
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
    gsl_vector_free(c);
    gsl_matrix_free(A);
    gsl_vector_free(b);
    var_arr_free(&var_arr);
    free(B);
    free(N);
    if (stream != stdin) {
        fclose(stream);
    }
    return 0;
}

void print_coefficient(double coef, uint32_t var_idx, int is_first) {
    char buf[PRINT_TERM_WIDTH + 1];

    if (fabs(coef) < 1e-9) {
        // Print spaces for alignment
        snprintf(buf, sizeof(buf), "%*s", PRINT_TERM_WIDTH, "");
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

    printf("%*s", PRINT_TERM_WIDTH, buf);
}

// Pretty print
void problem_print(const problem_t* problem_ptr, const char* name) {
    if (!problem_ptr || !name) {
        return;
    }

    printf("\n================== %s ==================\n", name);
    printf("%s z = ", problem_ptr->is_max ? "max" : "min");

    // Objective function
    const gsl_vector* c = problem_ptr->c;
    for (uint32_t i = 0; i < problem_ptr->m; i++) {
        double ci = gsl_vector_get(c, i);
        print_coefficient(ci, i, i == 0);
    }
    printf("\n\nconstraints:\n");

    // Constraints
    const gsl_matrix* A = problem_ptr->A;
    const gsl_vector* b = problem_ptr->b;
    for (uint32_t i = 0; i < problem_ptr->n; i++) {
        printf("\t");
        for (uint32_t j = 0; j < problem_ptr->m; j++) {
            double aij = gsl_matrix_get(A, i, j);
            print_coefficient(aij, j, j == 0);
        }

        double bi = gsl_vector_get(b, i);
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

uint32_t problem_is_milp(const problem_t* problem_ptr) {
    if (!problem_ptr) {
        return 0;
    }

    const var_arr_t* var_arr_ptr = &problem_ptr->var_arr;
    for (uint32_t i = 0; i < problem_ptr->m; i++) {
        const variable_t* v = var_arr_get(var_arr_ptr, i);
        if (variable_is_integer(v)) {
            return 1;
        }
    }

    return 0;
}

uint32_t problem_solve(problem_t* problem_ptr, solution_t* solution_ptr) {
    if (!problem_ptr || !solution_ptr) {
        fprintf(stderr, "Some arguments are NULL in problem_solve\n");
        return 0;
    }

    // Solve with B&B
    if (problem_is_milp(problem_ptr)) {
        return branch_and_bound(problem_ptr, solution_ptr);
    }

    // Solve with Primal Simplex
    uint32_t n = problem_ptr->n;
    uint32_t m = problem_ptr->m;
    uint32_t is_max = problem_ptr->is_max;

    gsl_vector_view c = gsl_vector_subvector(problem_ptr->c, 0, m);
    gsl_matrix_view A = gsl_matrix_submatrix(problem_ptr->A, 0, 0, n, m);
    gsl_vector_view b = gsl_vector_subvector(problem_ptr->b, 0, n);

    uint32_t iter_n = 0;
    uint32_t res = simplex_primal(n, m, is_max, &c.vector, &A.matrix, &b.vector, problem_ptr->B, problem_ptr->N,
                                  solution_ptr, &iter_n);

    solution_set_pI_iter(solution_ptr, problem_ptr->pI_iter);
    solution_set_pII_iter(solution_ptr, iter_n);

    return res;
}

void problem_free(problem_t* problem_ptr) {
    if (!problem_ptr) {
        return;
    }

    gsl_vector_free(problem_ptr->c);
    gsl_matrix_free(problem_ptr->A);
    gsl_vector_free(problem_ptr->b);
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

const gsl_vector* problem_c(const problem_t* problem_ptr) {
    return problem_ptr ? problem_ptr->c : NULL;
}

gsl_vector* problem_c_mut(problem_t* problem_ptr) {
    return problem_ptr ? problem_ptr->c : NULL;
}

const gsl_matrix* problem_A(const problem_t* problem_ptr) {
    return problem_ptr ? problem_ptr->A : NULL;
}

gsl_matrix* problem_A_mut(problem_t* problem_ptr) {
    return problem_ptr ? problem_ptr->A : NULL;
}

const gsl_vector* problem_b(const problem_t* problem_ptr) {
    return problem_ptr ? problem_ptr->b : NULL;
}

gsl_vector* problem_b_mut(problem_t* problem_ptr) {
    return problem_ptr ? problem_ptr->b : NULL;
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