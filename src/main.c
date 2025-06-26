#include "utils.h"
#include "simplex.h"
#include "branch_bound.h"
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_errno.h>

int read_model(FILE* stream, uint32_t* n, uint32_t* m, uint32_t* is_max, gsl_vector** c, gsl_matrix** A, gsl_vector** b,
               uint32_t** is_integer) {
    if (stream == stdin) {
        printf("n: ");
    }
    if (fscanf(stream, "%u", n) != 1 || *n > MAX_ROWS) {
        fprintf(stderr, "Invalid or too many rows (max %u).\n", MAX_ROWS);
        goto fail;
    }

    if (stream == stdin) {
        printf("m: ");
    }
    if (fscanf(stream, "%u", m) != 1) {
        fprintf(stderr, "Failed to read m.\n");
        goto fail;
    }
    fgetc(stream);  // consume newline

    if (stream == stdin) {
        printf("is_max: ");
    }
    if (fscanf(stream, "%u", is_max) != 1) {
        fprintf(stderr, "Failed to read is_max.\n");
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

    *is_integer = (uint32_t*)malloc(sizeof(uint32_t) * *m);
    if (!*is_integer) {
        fprintf(stderr, "Failed to allocate is_integer.\n");
        goto fail;
    }

    if (stream == stdin) {
        printf("is_integer (0/1 per variable): ");
    }

    for (uint32_t i = 0; i < *m; ++i) {
        if (fscanf(stream, "%u", &(*is_integer)[i]) != 1 || ((*is_integer)[i] != 0 && (*is_integer)[i] != 1)) {
            fprintf(stderr, "Invalid is_integer flag for variable %u.\n", i);
            goto fail;
        }
    }

    return 1;

fail:
    gsl_vector_free(*c);
    gsl_matrix_free(*A);
    gsl_vector_free(*b);
    free(*is_integer);
    return 0;
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
    uint32_t* is_integer = NULL;
    problem_t p = NULL;
    solution_t s = NULL;

    if (!read_model(stream, &n, &m, &is_max, &c, &A, &b, &is_integer)) {
        return EXIT_FAILURE;
    }

    if (read_from_file) {
        fclose(stream);
    }

    p = problem_new2(n, m, is_max, c, A, b, is_integer);
    if (!p) {
        fprintf(stderr, "Failed to create problem\n");
        gsl_vector_free(c);
        gsl_matrix_free(A);
        gsl_vector_free(b);
        free(is_integer);
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
    free(is_integer);
    problem_free(&p);
    solution_free(&s);

    return res;
}
