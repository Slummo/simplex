#include "solution.h"
#include <stdio.h>

struct solution {
    uint32_t n;             // Number of variables
    gsl_vector* x;          // Optimal solution (n)
    double z;               // Optimal value
    int32_t* basis;         // Final basis indices
    uint32_t is_unbounded;  // Boolean value to know if unbounded
    uint32_t pI_iter;       // Number of iterations of PhaseI to find a base
    uint32_t pII_iter;      // Number of iterations needed to solve
};

solution_t solution_new(uint32_t n, gsl_vector* x, int32_t* basis, uint32_t is_unbounded, uint32_t pI_iter,
                        uint32_t pII_iter) {
    solution_t s = (solution_t)malloc(sizeof(_solution));
    if (!s) {
        return NULL;
    }

    s->n = n;
    s->x = x;
    s->basis = basis;
    s->is_unbounded = is_unbounded;
    s->pI_iter = pI_iter;
    s->pII_iter = pII_iter;

    return s;
}

void solution_print(const solution_t s) {
    if (!s) {
        return;
    }

    printf("\n================== SOLUTION ==================\n");
    if (s->is_unbounded) {
        printf("infinite\n");
    } else {
        printf("Optimal found in %u iterations (PhaseI %u + PhaseII %u)\nz*: %lf\nx*: (", s->pI_iter + s->pII_iter,
               s->pI_iter, s->pII_iter, s->z);
    }

    if (!s->is_unbounded) {
        for (uint32_t i = 0; i < s->n; i++) {
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

gsl_vector* solution_optimal_vec(const solution_t s) {
    return s->x;
}

double solution_optimal_value(const solution_t s) {
    return s->z;
}

int32_t* solution_basis(const solution_t s) {
    return s->basis;
}

uint32_t is_solution_unbounded(const solution_t s) {
    return s->is_unbounded;
}

uint32_t solution_pI_iterations(const solution_t s) {
    return s->pI_iter;
}

uint32_t solution_pII_iterations(const solution_t s) {
    return s->pII_iter;
}

void solution_set_optimal_value(const solution_t s, double optimal_value) {
    s->z = optimal_value;
}