#include "problem.h"
#include "utils.h"
#include "simplex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/// Return the smaller of a and b
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

struct problem {
    uint32_t n;            // Number of constraints
    uint32_t m;            // Number of variables
    uint32_t is_max;       // Boolean value to know if its a maximization problem
    gsl_vector* c;         // Reduced costs (m)
    gsl_matrix* A;         // Constraints matrix (n x m)
    gsl_vector* b;         // RHS (n)
    int32_t* basis;        // Indices of basic variables (size n)
    int32_t* nonbasis;     // Indices of nonbasic variables (size m-n)
    uint32_t pI_iter;      // Number of iterations to find base with PhaseI
    uint32_t* is_integer;  // Boolean values to know if a variable is integer
};

// Duplicates each mallocable param
problem_t problem_new(uint32_t n, uint32_t m, uint32_t is_max, const gsl_vector* c, const gsl_matrix* A,
                      const gsl_vector* b, const int32_t* basis, uint32_t pI_iter, const uint32_t* is_integer) {
    if (!c || !A || !b || !basis) {
        return NULL;
    }

    problem_t p = (problem_t)malloc(sizeof(_problem));
    if (!p) {
        return NULL;
    }

    p->n = n;
    p->m = m;
    p->is_max = is_max;

    int32_t* used = NULL;

    // Duplicate vectors and matrix
    p->c = vector_duplicate(c);
    if (!p->c) {
        goto fail;
    }

    p->A = matrix_duplicate(A);
    if (!p->A) {
        goto fail;
    }

    p->b = vector_duplicate(b);
    if (!p->b) {
        goto fail;
    }

    // Copy basis array
    p->basis = (int32_t*)malloc(sizeof(int32_t) * n);
    if (!p->basis) {
        goto fail;
    }
    memcpy(p->basis, basis, sizeof(int32_t) * n);

    // Allocate and fill nonbasis array
    used = (int32_t*)calloc(m, sizeof(uint32_t));
    if (!used) {
        goto fail;
    }

    for (uint32_t i = 0; i < n; i++) {
        if (used[basis[i]]) {
            fprintf(stderr, "Duplicate basis index %u\n", basis[i]);
            goto fail;
        }
        used[basis[i]] = 1;
    }

    p->nonbasis = malloc(sizeof(int32_t) * (m - n));
    if (!p->nonbasis) {
        goto fail;
    }

    uint32_t nonbasis_count = 0;
    for (uint32_t i = 0; i < m; i++) {
        if (!used[i]) {
            p->nonbasis[nonbasis_count++] = i;
        }
    }
    free(used);
    used = NULL;

    p->pI_iter = pI_iter;

    // Handle integer flags if provided
    if (is_integer) {
        p->is_integer = (uint32_t*)malloc(sizeof(uint32_t) * m);
        if (!p->is_integer) {
            goto fail;
        }
        memcpy(p->is_integer, is_integer, sizeof(uint32_t) * m);
    } else {
        p->is_integer = NULL;
    }

    return p;

fail:
    gsl_vector_free(p->c);
    gsl_matrix_free(p->A);
    gsl_vector_free(p->b);
    free(p->basis);
    free(used);
    free(p->nonbasis);
    free(p->is_integer);
    problem_free(&p);
    return NULL;
}

problem_t problem_new2(uint32_t n, uint32_t m, uint32_t is_max, const gsl_vector* c, const gsl_matrix* A,
                       const gsl_vector* b, const uint32_t* is_integer) {
    // Find a base with phaseI
    uint32_t pI_iter = 0;
    int32_t* basis = simplex_phaseI(n, m, A, b, &pI_iter);
    if (!basis) {
        return NULL;
    }

    return problem_new(n, m, is_max, c, A, b, basis, pI_iter, is_integer);
}

problem_t problem_duplicate(const problem_t p) {
    if (!p) {
        return NULL;
    }

    return problem_new(p->n, p->m, p->is_max, p->c, p->A, p->b, p->basis, p->pI_iter, p->is_integer);
}

uint32_t problem_is_milp(const problem_t p) {
    if (!p) {
        return 0;
    }

    for (uint32_t i = 0; i < p->m; i++) {
        if (p->is_integer[i]) {
            return 1;
        }
    }

    return 0;
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

    printf("\n================== %s ==================\n", name);
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

uint32_t problem_pI_iter(const problem_t p) {
    return p ? p->pI_iter : 0;
}

const uint32_t* problem_integers(const problem_t p) {
    return p ? p->is_integer : NULL;
}

/* SETTERS */
void problem_set_basis_nonbasis(const problem_t p, int32_t* basis) {
    if (!p) {
        return;
    }

    p->basis = basis;

    p->nonbasis = (int32_t*)malloc(sizeof(int32_t) * (p->m - p->n));
    if (!p->nonbasis) {
        return;
    }

    // Fill non-basis from remaining indices
    uint32_t* used = calloc(p->m, sizeof(uint32_t));
    if (!used) {
        free(p->nonbasis);
        return;
    }

    for (uint32_t i = 0; i < p->n; i++) {
        if (used[basis[i]]) {
            fprintf(stderr, "Duplicate basis index %u\n", basis[i]);
            gsl_vector_free(p->c);
            gsl_matrix_free(p->A);
            gsl_vector_free(p->b);
            free(p->basis);
            free(p->nonbasis);
            free(p);
            free(used);
            return;
        }
        used[basis[i]] = 1;
    }

    // Fill non-basis with those not used
    uint32_t ni = 0;
    for (uint32_t i = 0; i < p->m; i++) {
        if (!used[i]) {
            p->nonbasis[ni++] = i;
        }
    }
    free(used);
}
