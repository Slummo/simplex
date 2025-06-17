#include <stdio.h>
#include <string.h>
#include "problem.h"
#include <utils.h>

int main(void) {
    size_t n, m;

    printf("n: ");
    scanf("%zu", &n);
    if (n > MAX_ROWS) {
        fprintf(stderr, "Too many rows (max %d).\n", MAX_ROWS);
        return EXIT_FAILURE;
    }

    printf("m: ");
    scanf("%zu", &m);
    getchar();  // consume leftover newline

    int is_max;
    printf("is_max: ");
    scanf("%d", &is_max);

    gsl_vector* c = read_vector(stdin, "c", m);
    if (!c) {
        return EXIT_FAILURE;
    }

    gsl_matrix* A = read_matrix(stdin, "A", n, m);
    if (!A) {
        gsl_vector_free(c);
        return EXIT_FAILURE;
    }

    gsl_vector* b = read_vector(stdin, "b", n);
    if (!b) {
        gsl_vector_free(c);
        gsl_matrix_free(A);
        return EXIT_FAILURE;
    }

    // Basis indices
    unsigned int* basis = (unsigned int*)malloc(sizeof(unsigned int) * n);
    if (!basis) {
        gsl_vector_free(c);
        gsl_matrix_free(A);
        gsl_vector_free(b);
        return EXIT_FAILURE;
    }

    printf("Enter indices of %zu basic variables (from 0 to %zu): ", n, m - 1);
    for (int i = 0; i < (int)n; i++) {
        // Dont check for < 0 since its unsigned and will get converted
        // to a number way over the value of m
        if (scanf("%u", &basis[i]) != 1 || basis[i] >= m) {
            fprintf(stderr, "Invalid basis index.\n");
            gsl_vector_free(c);
            gsl_matrix_free(A);
            gsl_vector_free(b);
            free(basis);
            return EXIT_FAILURE;
        }
    }
    getchar();

    problem_t p = problem_new(n, m, is_max, c, A, b, basis);
    if (!p) {
        fprintf(stderr, "Failed to create problem\n");
        gsl_vector_free(c);
        gsl_matrix_free(A);
        gsl_vector_free(b);
        free(basis);
        return EXIT_FAILURE;
    }

    solution_t s = solve(p);
    if (!s) {
        fprintf(stderr, "Failed to create solution\n");
        return EXIT_FAILURE;
    }

    solution_print(s);

    problem_free(&p);
    solution_free(&s);
    return EXIT_SUCCESS;
}
