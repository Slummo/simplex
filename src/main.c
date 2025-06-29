#include "rc.h"
#include "utils.h"
#include "simplex.h"
#include "branch_bound.h"
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_errno.h>

int read_model(FILE* stream, uint32_t* n, uint32_t* m, uint32_t* is_max, gsl_vector** c_raw, gsl_matrix** A_raw,
               gsl_vector** b_raw, variable_t*** variables_raw) {
    if (stream == stdin) {
        printf("n: ");
    }

    if (fscanf(stream, "%u", n) != 1 || *n > MAX_ROWS) {
        fprintf(stderr, "Invalid or too many rows (max %u)\n", MAX_ROWS);
        goto fail;
    }

    if (stream == stdin) {
        printf("m: ");
    }
    if (fscanf(stream, "%u", m) != 1) {
        fprintf(stderr, "Failed to read m\n");
        goto fail;
    }
    fgetc(stream);  // consume newline

    if (stream == stdin) {
        printf("is_max: ");
    }
    if (fscanf(stream, "%u", is_max) != 1) {
        fprintf(stderr, "Failed to read is_max\n");
        goto fail;
    }

    *c_raw = vector_read(stream, "c_raw", *m);
    if (!*c_raw) {
        goto fail;
    }

    *A_raw = matrix_read(stream, "A_raw", *n, *m);
    if (!*A_raw) {
        goto fail;
    }

    *b_raw = vector_read(stream, "b_raw", *n);
    if (!*b_raw) {
        goto fail;
    }

    *variables_raw = (variable_t**)malloc(sizeof(variable_t**) * *m);
    if (!*variables_raw) {
        fprintf(stderr, "Failed to allocate variables_raw\n");
        goto fail;
    }

    if (stream == stdin) {
        printf("type (0 = real / 1 = integer / 2 = binary) per variable: ");
    }

    uint32_t type = 0;
    for (uint32_t i = 0; i < *m; ++i) {
        if (fscanf(stream, "%u", &type) != 1 || (type != 0 && type != 1 && type != 2)) {
            fprintf(stderr, "Invalid type flag %u for variable %u\n", type, i);
            goto fail;
        } else {
            switch (type) {
                case 0: {
                    (*variables_raw)[i] = variable_new_real_positive(10e9);
                    break;
                }
                case 1: {
                    (*variables_raw)[i] = variable_new_integer_positive(10e9);
                    break;
                }
                case 2: {
                    (*variables_raw)[i] = variable_new_binary();
                    break;
                }
            }
        }
    }

    if (stream != stdin) {
        fclose(stream);
    }

    return 1;

fail:
    gsl_vector_free(*c_raw);
    gsl_matrix_free(*A_raw);
    gsl_vector_free(*b_raw);
    for (uint32_t i = 0; i < *m; i++) {
        variable_free(&(*variables_raw)[i]);
    }
    if (stream != stdin) {
        fclose(stream);
    }
    return 0;
}

int main(int argc, char** args) {
    if (argc > 2) {
        fprintf(stderr, "Wrong number of arguments! Usage: 'zmax <filename>' or just 'zmax' to read from stdin\n");
        return EXIT_FAILURE;
    }

    gsl_set_error_handler_off();

    FILE* stream = argc == 2 ? fopen(args[1], "r") : stdin;
    if (!stream) {
        perror("Failed to determine stream");
        return EXIT_FAILURE;
    }

    uint32_t n = 0;
    uint32_t m = 0;
    uint32_t is_max = 0;
    gsl_vector* c_raw = NULL;
    gsl_matrix* A_raw = NULL;
    gsl_vector* b_raw = NULL;
    variable_t** variables_raw = NULL;
    problem_t* p = NULL;
    solution_t* s = NULL;

    if (!read_model(stream, &n, &m, &is_max, &c_raw, &A_raw, &b_raw, &variables_raw)) {
        return EXIT_FAILURE;
    }

    p = problem_new2(n, m, is_max, c_raw, A_raw, b_raw, variables_raw);
    if (!p) {
        fprintf(stderr, "Failed to create problem\n");
        gsl_vector_free(c_raw);
        gsl_matrix_free(A_raw);
        gsl_vector_free(b_raw);
        for (uint32_t i = 0; i < m; i++) {
            variable_free(&variables_raw[i]);
        }
        return EXIT_FAILURE;
    }
    problem_print(p, "Problem");

    char* solver_name = NULL;

    if (problem_is_milp(p)) {
        solver_name = "Branch and bound";
        s = branch_and_bound(p);
    } else {
        solver_name = "Simplex";
        s = simplex_phaseII(p);
    }

    int res = EXIT_SUCCESS;
    if (!s) {
        fprintf(stderr, "Failed to solve with %s\n", solver_name);
        res = EXIT_FAILURE;
    }
    solution_print(s, "Solution");

    problem_free(&p);
    solution_free(&s);

    return res;
}
