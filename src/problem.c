#include "problem.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/// Return the smaller of a and b
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

struct problem {
    uint32_t n;         // Number of constraints
    uint32_t m;         // Number of variables
    uint32_t is_max;    // Boolean value to know if its a maximization problem
    gsl_vector* c;      // Reduced costs (m)
    gsl_matrix* A;      // Constraints matrix (n x m)
    gsl_vector* b;      // RHS (n)
    int32_t* basis;     // Indices of basic variables (size n)
    int32_t* nonbasis;  // Indices of nonbasic variables (size m-n)
};

// Frees arguments on error
problem_t problem_new(uint32_t n, uint32_t m, uint32_t is_max, gsl_vector* c, gsl_matrix* A, gsl_vector* b,
                      int32_t* basis) {
    if (!c || !A || !b || !basis) {
        return NULL;
    }

    problem_t p = (problem_t)malloc(sizeof(_problem));
    if (!p) {
        if (c)
            gsl_vector_free(c);
        if (A)
            gsl_matrix_free(A);
        if (b)
            gsl_vector_free(b);
        if (basis)
            free(basis);
        return NULL;
    }

    if (!is_max) {
        // Multiply by -1 to have the objective function
        // in standard form
        gsl_vector_scale(c, -1.0);
    }

    p->n = n;
    p->m = m;
    p->is_max = is_max;
    p->c = c;
    p->A = A;
    p->b = b;
    p->basis = basis;

    p->nonbasis = (int32_t*)malloc(sizeof(int32_t) * (m - n));
    if (!p->nonbasis) {
        problem_free(&p);
        return NULL;
    }

    // Fill non-basis from remaining indices
    uint32_t* used = calloc(m, sizeof(uint32_t));
    for (uint32_t i = 0; i < n; i++) {
        if (used[basis[i]]) {
            fprintf(stderr, "Duplicate basis index %u\n", basis[i]);
            problem_free(&p);
            return NULL;
        }
        used[basis[i]] = 1;
    }

    // Fill non-basis with those not used
    uint32_t ni = 0;
    for (uint32_t i = 0; i < m; i++) {
        if (!used[i]) {
            p->nonbasis[ni++] = i;
        }
    }
    free(used);

    return p;
}

problem_t problem_duplicate(const problem_t p) {
    if (!p) {
        return NULL;
    }

    gsl_vector* c = vector_duplicate(p->c);
    if (!c) {
        fprintf(stderr, "Error in problem_duplicate\n");
        return NULL;
    }

    gsl_matrix* A = matrix_duplicate(p->A);
    if (!A) {
        fprintf(stderr, "Error in problem_duplicate\n");
        gsl_vector_free(c);
        return NULL;
    }

    gsl_vector* b = vector_duplicate(p->b);
    if (!A) {
        fprintf(stderr, "Error in problem_duplicate\n");
        gsl_vector_free(c);
        gsl_matrix_free(A);
        return NULL;
    }

    int32_t* basis = (int32_t*)malloc(sizeof(int32_t) * p->n);
    if (!basis) {
        fprintf(stderr, "Error in problem_duplicate\n");
        gsl_vector_free(c);
        gsl_matrix_free(A);
        gsl_vector_free(b);
        return NULL;
    }

    for (uint32_t i = 0; i < p->n; i++) {
        basis[i] = p->basis[i];
    }

    return problem_new(p->n, p->m, p->is_max, c, A, b, basis);
}

#define TERM_WIDTH 8

void print_coefficient(double coef, uint32_t var_idx, int is_first) {
    char buf[TERM_WIDTH + 1];

    if (fabs(coef) < 1e-9) {
        // Print spaces for alignment
        snprintf(buf, sizeof(buf), "%*s", TERM_WIDTH, "");
    } else {
        char sign = coef < 0 ? '-' : (is_first ? ' ' : '+');
        double abs_val = fabs(coef);

        // Check if it's 1 or -1
        if (fabs(abs_val - 1.0) < 1e-9) {
            snprintf(buf, sizeof(buf), "%c x%-2u", sign, var_idx + 1);  // omit the 1
        } else if (fabs(abs_val - (int32_t)abs_val) < 1e-9) {
            snprintf(buf, sizeof(buf), "%c%2dx%-2u", sign, (int32_t)abs_val, var_idx + 1);
        } else {
            snprintf(buf, sizeof(buf), "%c%.1lfx%-2u", sign, abs_val, var_idx + 1);
        }
    }

    printf("%*s", TERM_WIDTH, buf);
}

// Pretty print
void problem_print(const problem_t p, const char* name) {
    if (!p || !name) {
        return;
    }

    printf("================== %s ==================\n", name);
    printf("%s z = ", p->is_max ? "max" : "min");

    // Objective function
    for (uint32_t i = 0; i < p->m; i++) {
        double ci = gsl_vector_get(p->c, i);
        print_coefficient(ci, i, i == 0);
    }
    printf("\n\nconstraints:\n");

    // Constraints
    for (uint32_t i = 0; i < p->n; i++) {
        printf("\t");
        for (uint32_t j = 0; j < p->m; j++) {
            double aij = gsl_matrix_get(p->A, i, j);
            print_coefficient(aij, j, j == 0);
        }

        double bi = gsl_vector_get(p->b, i);
        printf(" = ");
        if (fabs(bi - (int32_t)bi) < 1e-9) {
            printf("%d", (int32_t)bi);
        } else {
            printf("%.3lf", bi);
        }
        printf("\n");
    }

    printf("\n");
}

void problem_free(problem_t* pp) {
    if (!pp || !*pp) {
        return;
    }

    gsl_vector_free((*pp)->c);
    gsl_matrix_free((*pp)->A);
    gsl_vector_free((*pp)->b);
    free((*pp)->basis);
    free((*pp)->nonbasis);
    free(*pp);
    *pp = NULL;
}

/* GETTERS */

uint32_t problem_n(const problem_t p) {
    return p ? p->n : 0;
}

uint32_t problem_m(const problem_t p) {
    return p ? p->m : 0;
}

uint32_t problem_is_max(const problem_t p) {
    return p ? p->is_max : 0;
}

const gsl_vector* problem_c(const problem_t p) {
    return p ? p->c : NULL;
}

const gsl_matrix* problem_A(const problem_t p) {
    return p ? p->A : NULL;
}

const gsl_vector* problem_b(const problem_t p) {
    return p ? p->b : NULL;
}

const int32_t* problem_basis(const problem_t p) {
    return p ? p->basis : NULL;
}

int32_t* problem_basis_mut(const problem_t p) {
    return p ? p->basis : NULL;
}

const int32_t* problem_nonbasis(const problem_t p) {
    return p ? p->nonbasis : NULL;
}

int32_t* problem_nonbasis_mut(const problem_t p) {
    return p ? p->nonbasis : NULL;
}