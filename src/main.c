#include "utils.h"
#include "simplex.h"
#include "branch_bound.h"
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_errno.h>

int read_model(FILE* stream, uint32_t* n, uint32_t* m, uint32_t* is_max, gsl_vector** c, gsl_matrix** A,
               gsl_vector** b) {
    if (stream == stdin) {
        printf("n: ");
    }
    if (fscanf(stream, "%u", n) != 1 || *n > MAX_ROWS) {
        fprintf(stderr, "Invalid or too many rows (max %u).\n", MAX_ROWS);
        return 0;
    }

    if (stream == stdin) {
        printf("m: ");
    }
    if (fscanf(stream, "%u", m) != 1) {
        fprintf(stderr, "Failed to read m.\n");
        return 0;
    }
    fgetc(stream);  // consume newline

    if (stream == stdin) {
        printf("is_max: ");
    }
    if (fscanf(stream, "%u", is_max) != 1) {
        fprintf(stderr, "Failed to read is_max.\n");
        return 0;
    }

    *c = vector_read(stream, "c", *m);
    if (!*c) {
        return 0;
    }

    *A = matrix_read(stream, "A", *n, *m);
    if (!*A) {
        gsl_vector_free(*c);
        return 0;
    }

    *b = vector_read(stream, "b", *n);
    if (!*b) {
        gsl_vector_free(*c);
        gsl_matrix_free(*A);
        return 0;
    }

    return 1;
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

    uint32_t n, m, is_max;
    gsl_vector* c = NULL;
    gsl_matrix* A = NULL;
    gsl_vector* b = NULL;

    if (!read_model(stream, &n, &m, &is_max, &c, &A, &b)) {
        if (read_from_file) {
            fclose(stream);
        }

        gsl_vector_free(c);
        gsl_matrix_free(A);
        gsl_vector_free(b);

        return EXIT_FAILURE;
    }

    problem_t p;
    solution_t s;

    /* Simplex */
    p = problem_new2(n, m, is_max, c, A, b);
    if (!p) {
        fprintf(stderr, "Failed to create problem\n");
        return EXIT_FAILURE;
    }

    problem_print(p, "Problem");

    s = simplex_phaseII(p);
    if (!s) {
        fprintf(stderr, "Failed to solve with Simplex\n");
        problem_free(&p);
        return EXIT_FAILURE;
    }

    solution_print(s, "Simplex");
    solution_free(&s);
    problem_free(&p);

    /* Branch and bound */
    p = problem_new2(n, m, is_max, c, A, b);
    if (!p) {
        fprintf(stderr, "Failed to create problem\n");
        return EXIT_FAILURE;
    }

    problem_print(p, "Problem");

    s = branch_and_bound(p);
    if (!s) {
        fprintf(stderr, "Failed to solve with Branch and bound\n");
        problem_free(&p);
        return EXIT_FAILURE;
    }

    solution_print(s, "Branch and bound");
    solution_free(&s);
    problem_free(&p);

    gsl_vector_free(c);
    gsl_matrix_free(A);
    gsl_vector_free(b);
    return EXIT_SUCCESS;
}
