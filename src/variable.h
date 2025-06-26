#ifndef VARIABLE_H
#define VARIABLE_H

#include <stdint.h>

typedef enum { VAR_REAL, VAR_INTEGER, VAR_BINARY, VAR_ERR } variable_type_t;
const char* variable_type_to_str(variable_type_t vt);

typedef struct variable _variable, *variable_t;

variable_t variable_new(double lb, double ub, variable_type_t type);
variable_t variable_new_real_positive(double ub);
variable_t variable_new_integer_positive(double ub);
variable_t variable_new_binary();

void variable_print(const variable_t v);
void variable_free(variable_t* vp);
void variables_arr_free(variable_t** variables_arr_ptr, uint32_t var_num);

/* GETTERS */
double variable_lb(const variable_t v);
double variable_ub(const variable_t v);
uint32_t variable_is_integer(const variable_t v);

#endif