#ifndef VARIABLE_H
#define VARIABLE_H

#include <stdint.h>

/* VARIABLE */
typedef enum { VAR_REAL, VAR_INTEGER, VAR_BINARY, VAR_ERR } variable_type_t;
const char* variable_type_to_str(variable_type_t vt);

typedef struct variable variable_t;

variable_t* variable_new(double lb, double ub, variable_type_t type);
variable_t* variable_new_real_positive(double ub);
variable_t* variable_new_integer_positive(double ub);
variable_t* variable_new_binary();

variable_t* variable_duplicate(const variable_t* v);
void variable_print(const variable_t* v);
void variable_free(variable_t** vp);

/* Getters */
double variable_lb(const variable_t* v);
double variable_ub(const variable_t* v);
uint32_t variable_is_integer(const variable_t* v);

/* VARR */
// Wrapper for an array of variable_t* and its size
typedef struct varr varr_t;

/// @brief Creates a new varr_t struct. On fail the owner must free varr_raw
/// @param varr_raw A raw pointer to an array of variables
/// @param var_num The number of variables in the array
/// @return NULL on fail, a varr_t pointer on success
varr_t* varr_new(variable_t** varr_raw, uint32_t var_num);
const variable_t* varr_get(const varr_t* varr, uint32_t i);
const variable_t** varr_data(const varr_t* varr);
uint32_t varr_num(const varr_t* varr);
void varr_free(varr_t** varr_ptr);
void varr_drop(void* data);

#endif