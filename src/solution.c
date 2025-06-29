#include "solution.h"
#include "utils.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

struct solution {
    uint32_t n;             // Number of variables
    gsl_vector* x;          // Optimal solution (n)
    double z;               // Optimal value
    int32_t* basis;         // Final basis indices
    uint32_t is_unbounded;  // Boolean value to know if unbounded
    uint32_t pI_iter;       // Number of iterations of PhaseI to find a base
    uint32_t pII_iter;      // Number of iterations of PhaseII to find solution
};

solution_t* solution_new(uint32_t n, const gsl_vector* x, const int32_t* basis, uint32_t is_unbounded, uint32_t pI_iter,
                         uint32_t pII_iter) {
    if (!x || !basis) {
        return NULL;
    }

    solution_t* s = (solution_t*)malloc(sizeof(solution_t));
    if (!s) {
        return NULL;
    }

    s->n = n;
    s->x = vector_duplicate(x);
    if (!s->x) {
        free(s);
        return NULL;
    }

    s->basis = (int32_t*)malloc(sizeof(int32_t) * n);
    if (!s->basis) {
        gsl_vector_free(s->x);
        free(s);
        return NULL;
    }

    memcpy(s->basis, basis, sizeof(int32_t) * n);

    s->is_unbounded = is_unbounded;
    s->pI_iter = pI_iter;
    s->pII_iter = pII_iter;

    return s;
}

solution_t* solution_duplicate(const solution_t* s) {
    if (!s) {
        return NULL;
    }

    return solution_new(s->n, s->x, s->basis, s->is_unbounded, s->pI_iter, s->pII_iter);
}

// Checks if the i-th component of the solution is an integer
uint32_t solution_var_is_integer(const solution_t* s, uint32_t i) {
    if (!s) {
        fprintf(stderr, "s is NULL in solution_var_is_integer\n");
        return 0;
    }

    double xi = gsl_vector_get(s->x, i);
    double diff = fabs(xi - round(xi));
    return diff < 1e-8;
}

// Pretty print
void solution_print(const solution_t* s, const char* name) {
    if (!s) {
        return;
    }

    printf("\n================== %s ==================\n", name);
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

void solution_free(solution_t** sp) {
    if (!sp || !*sp) {
        return;
    }

    gsl_vector_free((*sp)->x);
    free((*sp)->basis);
    free(*sp);
    *sp = NULL;
}

/* GETTERS */

const gsl_vector* solution_x(const solution_t* s) {
    return s ? s->x : NULL;
}

gsl_vector* solution_x_mut(const solution_t* s) {
    return s ? s->x : NULL;
}

double solution_z(const solution_t* s) {
    return s ? s->z : 0.0;
}

const int32_t* solution_basis(const solution_t* s) {
    return s ? s->basis : NULL;
}

uint32_t solution_is_unbounded(const solution_t* s) {
    return s ? s->is_unbounded : 0;
}

uint32_t solution_pI_iterations(const solution_t* s) {
    return s ? s->pI_iter : 0;
}

uint32_t solution_pII_iterations(const solution_t* s) {
    return s ? s->pII_iter : 0;
}

/* SETTERS */

void solution_set_optimal_value(solution_t* s, double optimal_value) {
    if (!s) {
        return;
    }

    s->z = optimal_value;
}