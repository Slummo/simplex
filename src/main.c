#include "rc.h"
#include "utils.h"
#include "simplex.h"
#include "branch_bound.h"
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_errno.h>

int read_model(FILE* stream, uint32_t* n, uint32_t* m, uint32_t* is_max, gsl_vector** c, gsl_matrix** A, gsl_vector** b,
               variable_t** variables) {
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

    *c = vector_read(stream, "c", *m);
    if (!*c) {
        goto fail;
    }

    *A = matrix_read(stream, "A", *n, *m);
    if (!*A) {
        goto fail;
    }

    *b = vector_read(stream, "b", *n);
    if (!*b) {
        goto fail;
    }

    *variables = (variable_t*)malloc(sizeof(variable_t) * *m);
    if (!*variables) {
        fprintf(stderr, "Failed to allocate variables\n");
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
                    (*variables)[i] = variable_new_real_positive(10e9);
                    break;
                }
                case 1: {
                    (*variables)[i] = variable_new_integer_positive(10e9);
                    break;
                }
                case 2: {
                    (*variables)[i] = variable_new_binary();
                    break;
                }
            }
        }
    }

    return 1;

fail:
    gsl_vector_free(*c);
    gsl_matrix_free(*A);
    gsl_vector_free(*b);
    variables_arr_free(variables, *m);
    return 0;
}

void drop_int_ptr(void* ptr) {
    free(ptr);
}

int main(int argc, char** args) {
    if (argc > 2) {
        fprintf(stderr, "Wrong number of arguments! Usage: 'zmax <filename>' or just 'zmax' to read from stdin\n");
        return EXIT_FAILURE;
    }

    FILE* stream = NULL;
    uint32_t read_from_file = argc == 2;
    if (read_from_file) {
        stream = fopen(args[1], "r");
        if (!stream) {
            perror("Failed to open file");
            return EXIT_FAILURE;
        }
    } else {
        stream = stdin;
    }

    gsl_set_error_handler_off();

    uint32_t n = 0;
    uint32_t m = 0;
    uint32_t is_max = 0;
    gsl_vector* c = NULL;
    gsl_matrix* A = NULL;
    gsl_vector* b = NULL;
    variable_t* variables = NULL;
    problem_t p = NULL;
    solution_t s = NULL;

    if (!read_model(stream, &n, &m, &is_max, &c, &A, &b, &variables)) {
        return EXIT_FAILURE;
    }

    if (read_from_file) {
        fclose(stream);
    }

    p = problem_new2(n, m, is_max, c, A, b, variables);
    if (!p) {
        fprintf(stderr, "Failed to create problem\n");
        gsl_vector_free(c);
        gsl_matrix_free(A);
        gsl_vector_free(b);
        variables_arr_free(&variables, m);
        return EXIT_FAILURE;
    }
    problem_print(p, "Problem");

    char* solver_name = NULL;
    int res = EXIT_SUCCESS;

    if (problem_is_milp(p)) {
        solver_name = "Branch and bound";
        s = branch_and_bound(p);
    } else {
        solver_name = "Simplex";
        s = simplex_phaseII(p);
    }

    if (!s) {
        fprintf(stderr, "Failed to solve with %s\n", solver_name);
        res = EXIT_FAILURE;
    }
    solution_print(s, "Solution");

    gsl_vector_free(c);
    gsl_matrix_free(A);
    gsl_vector_free(b);
    variables_arr_free(&variables, m);
    problem_free(&p);
    solution_free(&s);

    return res;
}
