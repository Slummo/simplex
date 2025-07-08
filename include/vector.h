#ifndef VECTOR_H
#define VECTOR_H

#include <stdint.h>
#include <gsl/gsl_vector.h>

typedef struct vector {
    double* data;
    uint32_t capacity;
} vector_t;

uint32_t vector_init(vector_t* vector_ptr, uint32_t capacity);
uint32_t vector_from_stream(vector_t* vector_ptr, FILE* stream, char* name, uint32_t capacity, uint32_t size);
uint32_t vector_as_gsl_view(vector_t* vector_ptr, uint32_t offset, uint32_t length, gsl_vector_view* view_ptr);
double vector_get(const vector_t* vector_ptr, uint32_t i);
uint32_t vector_set(vector_t* vector_ptr, uint32_t i, double value);
void vector_free(vector_t* vector_ptr);

#endif