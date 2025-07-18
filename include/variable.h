#ifndef VARIABLE_H
#define VARIABLE_H

#include <stdio.h>
#include <stdint.h>

/* VARIABLE */
typedef enum { VAR_REAL, VAR_INTEGER, VAR_BINARY, VAR_ERR } variable_type_t;
const char* variable_type_to_str(variable_type_t vt);

typedef struct variable {
    double lb;
    double ub;
    variable_type_t type;
} variable_t;

uint32_t variable_init(variable_t* variable_ptr, double lb, double ub, variable_type_t type);
uint32_t variable_init_real_positive(variable_t* variable_ptr, double ub);
uint32_t variable_init_integer_positive(variable_t* variable_ptr, double ub);
uint32_t variable_init_binary(variable_t* variable_ptr);

variable_t variable_copy(variable_t other);

uint32_t variable_is_real(const variable_t* variable_ptr);
uint32_t variable_is_integer(const variable_t* variable_ptr);
uint32_t variable_is_binary(const variable_t* variable_ptr);

void variable_print(const variable_t* v);
void variable_free(variable_t* variable_ptr);

/* VAR_ARRAY */
typedef struct var_arr {
    variable_t* data;
    uint32_t length;
    uint32_t capacity;
} var_arr_t;

uint32_t var_arr_init(var_arr_t* var_arr_ptr, uint32_t initial_capacity);
uint32_t var_arr_from_stream(var_arr_t* var_arr_ptr, FILE* stream, uint32_t capacity, uint32_t size);
const variable_t* var_arr_get(const var_arr_t* array_ptr, uint32_t i);
uint32_t var_arr_length(const var_arr_t* array_ptr);
uint32_t var_arr_capacity(const var_arr_t* array_ptr);
uint32_t var_arr_push(var_arr_t* array_ptr, variable_t* variable_ptr);
uint32_t var_arr_duplicate(const var_arr_t* original, var_arr_t* new);
void var_arr_free(var_arr_t* array_ptr);

#endif