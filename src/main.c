#include <stdio.h>
#include <string.h>
#include "problem.h"
#include <utils.h>

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
        fprintf(stderr,
                "Wrong number of arguments! Usage: 'simplex <filename>' or just 'simplex' to read from stdin\n");
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

    uint32_t p1_iter = 0;
    uint32_t* basis = problem_find_basis(n, m, A, b, &p1_iter);
    if (!basis) {
        gsl_vector_free(c);
        gsl_matrix_free(A);
        gsl_vector_free(b);
        return EXIT_FAILURE;
    }

    problem_t p = problem_new(n, m, is_max, c, A, b, basis);
    if (!p) {
        fprintf(stderr, "Failed to create problem\n");
        return EXIT_FAILURE;
    }

    problem_print(p, "Problem");

    solution_t s = solve(p);
    if (!s) {
        fprintf(stderr, "Failed to create solution\n");
        problem_free(&p);
        return EXIT_FAILURE;
    }
    // Update with iterations from Phase1
    solution_set_iterations(s, p1_iter + solution_iterations(s));

    solution_print(s);

    problem_free(&p);
    solution_free(&s);
    return EXIT_SUCCESS;
}
