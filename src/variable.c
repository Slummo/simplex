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

uint32_t variable_init(variable_t* variable_ptr, double lb, double ub, variable_type_t type) {
    if (!variable_ptr) {
        return 0;
    }

    variable_ptr->lb = lb;
    variable_ptr->ub = ub;
    variable_ptr->type = type;

    return 1;
}

uint32_t variable_init_real_positive(variable_t* variable_ptr, double ub) {
    return variable_init(variable_ptr, 0.0, ub, VAR_REAL);
}

uint32_t variable_init_integer_positive(variable_t* variable_ptr, double ub) {
    return variable_init(variable_ptr, 0.0, ub, VAR_INTEGER);
}

uint32_t variable_init_binary(variable_t* variable_ptr) {
    return variable_init(variable_ptr, 0.0, 1.0, VAR_BINARY);
}

variable_t variable_copy(variable_t other) {
    return other;
}

void variable_print(const variable_t* v) {
    if (!v) {
        return;
    }

    printf("Variable[lb: %.3lf, ub: %.3lf, type: %s]\n", v->lb, v->ub, variable_type_to_str(v->type));
}

void variable_free(variable_t* variable_ptr) {
    if (!variable_ptr) {
        return;
    }

    return;
}

/* VARR */
uint32_t var_arr_init(var_arr_t* var_arr_ptr, uint32_t initial_capacity) {
    if (!var_arr_ptr) {
        return 0;
    }

    var_arr_ptr->data = (variable_t*)malloc(sizeof(variable_t) * initial_capacity);
    if (!var_arr_ptr->data) {
        return 0;
    }
    var_arr_ptr->length = 0;
    var_arr_ptr->capacity = initial_capacity;

    return 1;
}

uint32_t var_arr_from_stream(var_arr_t* var_arr_ptr, FILE* stream, uint32_t capacity, uint32_t size) {
    if (size > capacity) {
        fprintf(stderr, "Requested size exceeds capacity in var_arr_from_stream\n");
        return 0;
    }

    if (!var_arr_init(var_arr_ptr, capacity)) {
        fprintf(stderr, "Failed to init var_arr in var_arr_from_stream\n");
        return 0;
    }

    if (stream == stdin) {
        printf("type (0 = real / 1 = integer / 2 = binary) per variable: ");
        fflush(stdout);
    }

    uint32_t type = 0;
    variable_t v;
    for (uint32_t i = 0; i < size; ++i) {
        if (fscanf(stream, "%u", &type) != 1 || (type != 0 && type != 1 && type != 2)) {
            fprintf(stderr, "Invalid type flag %u for variable %u\n", type, i);
            goto fail;
        } else {
            switch (type) {
                case 0: {
                    if (!variable_init_real_positive(&v, 10e9) || !var_arr_push(var_arr_ptr, &v)) {
                        goto fail;
                    }
                    break;
                }
                case 1: {
                    if (!variable_init_integer_positive(&v, 10e9) || !var_arr_push(var_arr_ptr, &v)) {
                        goto fail;
                    }
                    break;
                }
                case 2: {
                    if (!variable_init_binary(&v) || !var_arr_push(var_arr_ptr, &v)) {
                        goto fail;
                    }
                    break;
                }
            }
        }
    }

    return 1;

fail:
    var_arr_free(var_arr_ptr);
    return 0;
}

const variable_t* var_arr_get(const var_arr_t* array_ptr, uint32_t i) {
    if (array_ptr && i < array_ptr->length && array_ptr->data) {
        return &array_ptr->data[i];
    }

    __builtin_trap();
    return NULL;
}

// Takes ownership of *variable_ptr
uint32_t var_arr_push(var_arr_t* array_ptr, variable_t* variable_ptr) {
    if (!array_ptr || array_ptr->length >= array_ptr->capacity || !variable_ptr) {
        return 0;
    }

    array_ptr->data[array_ptr->length++] = variable_copy(*variable_ptr);
    return 1;
}

void var_arr_free(var_arr_t* array_ptr) {
    if (!array_ptr || !array_ptr->data) {
        return;
    }

    for (uint32_t i = 0; i < array_ptr->length; i++) {
        variable_free((variable_t*)var_arr_get(array_ptr, i));
    }

    free(array_ptr->data);
    array_ptr->data = NULL;
    array_ptr->length = 0;
    array_ptr->capacity = 0;
}