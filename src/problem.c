#include "problem.h"
#include "utils.h"
#include "simplex.h"
#include "rc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/// Return the smaller of a and b
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

struct problem {
    uint32_t n;        // Number of constraints
    uint32_t m;        // Number of variables
    uint32_t is_max;   // Boolean value to know if its a maximization problem
    rc_t c;            // Reduced costs (m)
    rc_t A;            // Constraints matrix (n x m)
    rc_t b;            // RHS (n)
    rc_t basis;        // Indices of basic variables (size n)
    rc_t nonbasis;     // Indices of nonbasic variables (size m-n)
    uint32_t pI_iter;  // Number of iterations to find base with PhaseI
    rc_t varr;         // Array of variables
};

void gsl_vector_drop(void* data) {
    gsl_vector_free((gsl_vector*)data);
}

void gsl_matrix_drop(void* data) {
    gsl_matrix_free((gsl_matrix*)data);
}

problem_t* problem_new(uint32_t n, uint32_t m, uint32_t is_max, const gsl_vector* c_raw, const gsl_matrix* A_raw,
                       const gsl_vector* b_raw, const int32_t* basis_raw, uint32_t pI_iter,
                       variable_t** variables_raw) {
    if (!c_raw || !A_raw || !b_raw || !basis_raw || !variables_raw) {
        return NULL;
    }

    problem_t* p = (problem_t*)malloc(sizeof(problem_t));
    if (!p) {
        return NULL;
    }

    p->n = n;
    p->m = m;
    p->is_max = is_max;

    int32_t* used = NULL;
    int32_t* nonbasis = NULL;
    varr_t* varr = NULL;

    p->c = rc_new((void*)c_raw, gsl_vector_drop);
    if (!p->c) {
        goto fail;
    }

    p->A = rc_new((void*)A_raw, gsl_matrix_drop);
    if (!p->A) {
        goto fail;
    }

    p->b = rc_new((void*)b_raw, gsl_vector_drop);
    if (!p->b) {
        goto fail;
    }

    p->basis = rc_new((void*)basis_raw, free);
    if (!p->basis) {
        goto fail;
    }

    used = (int32_t*)calloc(m, sizeof(uint32_t));
    if (!used) {
        goto fail;
    }

    for (uint32_t i = 0; i < n; i++) {
        if (used[basis_raw[i]]) {
            fprintf(stderr, "Duplicate basis index %u\n", basis_raw[i]);
            goto fail;
        }
        used[basis_raw[i]] = 1;
    }

    nonbasis = (int32_t*)malloc(sizeof(int32_t) * (m - n));
    if (!nonbasis) {
        goto fail;
    }

    uint32_t nonbasis_count = 0;
    for (uint32_t i = 0; i < m; i++) {
        if (!used[i]) {
            nonbasis[nonbasis_count++] = i;
        }
    }
    free(used);
    used = NULL;

    p->nonbasis = rc_new((void*)nonbasis, free);
    if (!p->nonbasis) {
        goto fail;
    }

    p->pI_iter = pI_iter;

    varr = varr_new(variables_raw, m);
    if (!varr) {
        goto fail;
    }

    p->varr = rc_new((void*)varr, varr_drop);
    if (!p->varr) {
        goto fail;
    }

    return p;

fail:
    rc_free(&p->c);
    rc_free(&p->A);
    rc_free(&p->b);
    rc_free(&p->basis);
    free(used);
    free(nonbasis);
    rc_free(&p->nonbasis);
    varr_free(&varr);
    rc_free(&p->varr);
    problem_free(&p);
    return NULL;
}

problem_t* problem_new2(uint32_t n, uint32_t m, uint32_t is_max, const gsl_vector* c_raw, const gsl_matrix* A_raw,
                        const gsl_vector* b_raw, variable_t** variables_raw) {
    // Find a base with phaseI
    uint32_t pI_iter = 0;
    int32_t* basis = simplex_phaseI(n, m, A_raw, b_raw, variables_raw, &pI_iter);
    if (!basis) {
        return NULL;
    }

    return problem_new(n, m, is_max, c_raw, A_raw, b_raw, basis, pI_iter, variables_raw);
}

problem_t* problem_duplicate(const problem_t* p) {
    if (!p || !p->c || !p->A || !p->b || !p->basis || !p->nonbasis || !p->varr) {
        return NULL;
    }

    problem_t* p2 = (problem_t*)malloc(sizeof(problem_t));
    if (!p2) {
        return NULL;
    }

    // Shallow copy
    *p2 = *p;

    // Clone
    p2->c = rc_clone(p->c);
    p2->A = rc_clone(p->A);
    p2->b = rc_clone(p->b);
    p2->basis = rc_clone(p->basis);
    p2->nonbasis = rc_clone(p->nonbasis);
    p2->varr = rc_clone(p->varr);

    return p2;
}

uint32_t problem_is_milp(const problem_t* p) {
    if (!p) {
        return 0;
    }

    const varr_t* varr = problem_varr(p);
    for (uint32_t i = 0; i < p->m; i++) {
        if (variable_is_integer(varr_get(varr, i))) {
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
void problem_print(const problem_t* p, const char* name) {
    if (!p || !name) {
        return;
    }

    printf("\n================== %s ==================\n", name);
    printf("%s z = ", p->is_max ? "max" : "min");

    // Objective function
    const gsl_vector* c = problem_c(p);
    for (uint32_t i = 0; i < p->m; i++) {
        double ci = gsl_vector_get(c, i);
        print_coefficient(ci, i, i == 0);
    }
    printf("\n\nconstraints:\n");

    // Constraints
    const gsl_matrix* A = problem_A(p);
    const gsl_vector* b = problem_b(p);
    for (uint32_t i = 0; i < p->n; i++) {
        printf("\t");
        for (uint32_t j = 0; j < p->m; j++) {
            double aij = gsl_matrix_get(A, i, j);
            print_coefficient(aij, j, j == 0);
        }

        double bi = gsl_vector_get(b, i);
        printf(" = ");
        if (fabs(bi - (int32_t)bi) < 1e-9) {
            printf("%d", (int32_t)bi);
        } else {
            printf("%.3lf", bi);
        }
        printf("\n");
    }

    printf("\nvariables:\n");
    const varr_t* varr = problem_varr(p);
    for (uint32_t i = 0; i < p->m; i++) {
        printf("\t");
        variable_print(varr_get(varr, i));
    }
    puts("");
}

void problem_free(problem_t** pp) {
    if (!pp || !*pp) {
        return;
    }

    rc_free(&(*pp)->c);
    rc_free(&(*pp)->A);
    rc_free(&(*pp)->b);
    rc_free(&(*pp)->basis);
    rc_free(&(*pp)->nonbasis);
    rc_free(&(*pp)->varr);
    free(*pp);
    *pp = NULL;
}

/* GETTERS */

uint32_t problem_n(const problem_t* p) {
    return p ? p->n : 0;
}

uint32_t problem_m(const problem_t* p) {
    return p ? p->m : 0;
}

uint32_t problem_is_max(const problem_t* p) {
    return p ? p->is_max : 0;
}

const gsl_vector* problem_c(const problem_t* p) {
    return p ? (const gsl_vector*)rc_data(p->c) : NULL;
}

const gsl_matrix* problem_A(const problem_t* p) {
    return p ? (const gsl_matrix*)rc_data(p->A) : NULL;
}

const gsl_vector* problem_b(const problem_t* p) {
    return p ? (const gsl_vector*)rc_data(p->b) : NULL;
}

const int32_t* problem_basis(const problem_t* p) {
    return p ? (const int32_t*)rc_data(p->basis) : NULL;
}

int32_t* problem_basis_mut(const problem_t* p) {
    return p ? (int32_t*)rc_data_mut(p->basis) : NULL;
}

const int32_t* problem_nonbasis(const problem_t* p) {
    return p ? (const int32_t*)rc_data(p->nonbasis) : NULL;
}

int32_t* problem_nonbasis_mut(const problem_t* p) {
    return p ? (int32_t*)rc_data_mut(p->nonbasis) : NULL;
}

uint32_t problem_pI_iter(const problem_t* p) {
    return p ? p->pI_iter : 0;
}

const varr_t* problem_varr(const problem_t* p) {
    return p ? (varr_t*)rc_data(p->varr) : NULL;
}