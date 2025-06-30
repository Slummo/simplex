#include "variable.h"
#include <stdio.h>
#include <stdlib.h>

/* VARIABLE */

const char* variable_type_to_str(variable_type_t vt) {
    switch (vt) {
        case VAR_REAL: {
            return "VAR_REAL";
        }
        case VAR_INTEGER: {
            return "VAR_INTEGER";
        }
        case VAR_BINARY: {
            return "VAR_BINARY";
        }
        case VAR_ERR: {
            return "VAR_ERR";
        }
        default: {
            return "VAR_UNKNOWN";
        }
    }
}

struct variable {
    double lb;
    double ub;
    variable_type_t type;
};

variable_t* variable_new(double lb, double ub, variable_type_t type) {
    variable_t* v = (variable_t*)malloc(sizeof(variable_t));
    if (!v) {
        return NULL;
    }

    v->lb = lb;
    v->ub = ub;
    v->type = type;

    return v;
}

variable_t* variable_new_real_positive(double ub) {
    return variable_new(0, ub, VAR_REAL);
}

variable_t* variable_new_integer_positive(double ub) {
    return variable_new(0, ub, VAR_INTEGER);
}

variable_t* variable_new_binary() {
    return variable_new(0, 1, VAR_BINARY);
}

variable_t* variable_duplicate(const variable_t* v) {
    return v ? variable_new(v->lb, v->ub, v->type) : NULL;
}

void variable_print(const variable_t* v) {
    if (!v) {
        return;
    }

    printf("Variable[lb: %.3lf, ub: %.3lf, type: %s]\n", v->lb, v->ub, variable_type_to_str(v->type));
}

void variable_free(variable_t** vp) {
    if (!vp || !*vp) {
        return;
    }

    free(*vp);
    *vp = NULL;
}

/* Getters */

double variable_lb(const variable_t* v) {
    return v ? v->lb : -1.0;
}

double variable_ub(const variable_t* v) {
    return v ? v->ub : -1.0;
}

uint32_t variable_is_integer(const variable_t* v) {
    return v ? v->type == VAR_INTEGER : 0;
}

/* VARR */

struct varr {
    variable_t** data;
    uint32_t n;
};

varr_t* varr_new(variable_t** varr_raw, uint32_t var_num) {
    if (!varr_raw || !var_num) {
        return NULL;
    }

    varr_t* varr = (varr_t*)malloc(sizeof(varr_t));
    if (!varr) {
        return NULL;
    }

    varr->data = varr_raw;
    varr->n = var_num;

    return varr;
}

const variable_t* varr_get(const varr_t* varr, uint32_t i) {
    return varr ? varr->data[i] : NULL;
}

const variable_t** varr_data(const varr_t* varr) {
    return varr ? (const variable_t**)varr->data : NULL;
}

uint32_t varr_num(const varr_t* varr) {
    return varr ? varr->n : 0;
}

void varr_free(varr_t** varr_ptr) {
    if (!varr_ptr || !*varr_ptr) {
        return;
    }

    for (uint32_t i = 0; i < (*varr_ptr)->n; i++) {
        variable_free(&(*varr_ptr)->data[i]);
    }

    free((*varr_ptr)->data);
    free(*varr_ptr);
    *varr_ptr = NULL;
}

void varr_drop(void* data) {
    if (!data) {
        return;
    }

    varr_t* varr = (varr_t*)data;

    for (uint32_t i = 0; i < varr->n; i++) {
        variable_free(&varr->data[i]);
    }

    free(varr->data);
    free(varr);
}