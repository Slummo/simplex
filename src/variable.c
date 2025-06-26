#include "variable.h"
#include <stdio.h>
#include <stdlib.h>

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

variable_t variable_new(double lb, double ub, variable_type_t type) {
    variable_t v = (variable_t)malloc(sizeof(_variable));
    if (!v) {
        return NULL;
    }

    v->lb = lb;
    v->ub = ub;
    v->type = type;

    return v;
}

variable_t variable_new_real_positive(double ub) {
    return variable_new(0, ub, VAR_REAL);
}

variable_t variable_new_integer_positive(double ub) {
    return variable_new(0, ub, VAR_INTEGER);
}

variable_t variable_new_binary() {
    return variable_new(0, 1, VAR_BINARY);
}

void variable_print(const variable_t v) {
    if (!v) {
        return;
    }

    printf("Variable[lb: %.3lf, ub: %.3lf, type: %s]\n", v->lb, v->ub, variable_type_to_str(v->type));
}

void variable_free(variable_t* vp) {
    if (*vp || !*vp) {
        return;
    }

    free(*vp);
    *vp = NULL;
}

void variables_arr_free(variable_t** variables_arr_ptr, uint32_t var_num) {
    if (!variables_arr_ptr || !*variables_arr_ptr) {
        return;
    }

    for (uint32_t i = 0; i < var_num; i++) {
        variable_free(&(*variables_arr_ptr)[i]);
    }
    free(*variables_arr_ptr);
    *variables_arr_ptr = NULL;
}

/* GETTERS */

double variable_lb(const variable_t v) {
    return v ? v->lb : -1.0;
}

double variable_ub(const variable_t v) {
    return v ? v->ub : -1.0;
}

uint32_t variable_is_integer(const variable_t v) {
    return v ? v->type == VAR_INTEGER : 0;
}