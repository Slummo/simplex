#include <stdio.h>
#include <string.h>
#include "problem.h"

int main(void) {
    int n, m;

    printf("n: ");
    scanf("%d", &n);
    if (n > MAX_ROWS) {
        fprintf(stderr, "Too many rows (max %d).\n", MAX_ROWS);
        return EXIT_FAILURE;
    }

    printf("m: ");
    scanf("%d", &m);
    getchar();  // consume leftover newline

    // Vector c
    gsl_vector* c = gsl_vector_alloc(m);
    if (!c) {
        fprintf(stderr, "Failed to allocate vector c\n");
        return EXIT_FAILURE;
    }

    printf("Enter %d values for vector c: ", m);
    if (gsl_vector_fscanf(stdin, c) != 0) {
        fprintf(stderr, "Failed to read vector c\n");
        gsl_vector_free(c);
        return EXIT_FAILURE;
    }

    // Matrix A
    gsl_matrix* A = gsl_matrix_alloc(n, m);
    if (!A) {
        fprintf(stderr, "Failed to allocate matrix A\n");
        gsl_vector_free(c);
        return EXIT_FAILURE;
    }

    printf("Enter %dx%d = %d values for matrix A (row-major): ", n, m, n * m);
    if (gsl_matrix_fscanf(stdin, A) != 0) {
        fprintf(stderr, "Failed to read matrix A\n");
        gsl_vector_free(c);
        gsl_matrix_free(A);
        return EXIT_FAILURE;
    }

    // Vector b
    gsl_vector* b = gsl_vector_alloc(n);
    if (!b) {
        fprintf(stderr, "Failed to allocate vector b\n");
        gsl_vector_free(c);
        gsl_matrix_free(A);
        return EXIT_FAILURE;
    }

    printf("Enter %d values for vector b: ", n);
    if (gsl_vector_fscanf(stdin, b) != 0) {
        fprintf(stderr, "Failed to read vector b\n");
        gsl_vector_free(c);
        gsl_matrix_free(A);
        gsl_vector_free(b);
        return EXIT_FAILURE;
    }

    // Basis indices
    int* basis = (int*)malloc(sizeof(int) * n);
    printf("Enter indices of %d basic variables (space-separated, from 0 to %d): ", n, m - 1);
    for (int i = 0; i < n; i++) {
        if (scanf("%d", &basis[i]) != 1 || basis[i] < 0 || basis[i] >= m) {
            fprintf(stderr, "Invalid basis index.\n");
            return EXIT_FAILURE;
        }
    }
    getchar();

    problem_t p = problem_new(n, m, c, A, b, basis);

    solution_t s = solve(p);
    solution_print(s);

    problem_free(&p);
    solution_free(&s);
    return EXIT_SUCCESS;
}
