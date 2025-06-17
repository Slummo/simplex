#include "solution.h"
#include <stdio.h>

struct solution {
    int n;             // Number of variables
    gsl_vector* x;     // Optimal solution (n)
    double z;          // Optimal value
    int is_unbounded;  // Boolean value to know if unbounded
    int n_iter;        // Number of iterations needed to solve
};

solution_t solution_new(int n, gsl_vector* x, double z, int is_unbounded, int n_iter) {
    solution_t s = (solution_t)malloc(sizeof(_solution));
    if (!s) {
        return NULL;
    }

    s->n = n;
    s->x = x;
    s->z = z;
    s->is_unbounded = is_unbounded;
    s->n_iter = n_iter;

    return s;
}

void solution_print(solution_t s) {
    if (!s) {
        return;
    }

    printf("\n==================SOLUTION==================\n");
    if (s->is_unbounded) {
        printf("infinite\n");
    } else {
        printf("Optimal found in %d iterations\nz*: %lf\nx*: (", s->n_iter, s->z);
    }

    if (!s->is_unbounded) {
        for (int i = 0; i < s->n; i++) {
            printf("%.3lf", gsl_vector_get(s->x, i));
            if (i < s->n - 1) {
                printf(", ");
            }
        }
        puts(")");
    }
}

void solution_free(solution_t* sp) {
    if (!sp || !*sp) {
        return;
    }

    gsl_vector_free((*sp)->x);
    free(*sp);
    *sp = NULL;
}

gsl_vector* solution_optimal_vec_ptr(const solution_t s) {
    return s->x;
}

double* solution_optimal_value_ptr(const solution_t s) {
    return &s->z;
}

int is_solution_unbounded(const solution_t s) {
    return s->is_unbounded;
}