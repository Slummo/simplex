#include "vector.h"
#include <stdio.h>
#include <stdlib.h>

uint32_t vector_init(vector_t* vector_ptr, uint32_t capacity) {
    if (!vector_ptr) {
        return 0;
    }

    vector_ptr->data = (double*)calloc(capacity, sizeof(double));
    if (!vector_ptr->data) {
        return 0;
    }

    vector_ptr->capacity = capacity;
    return 1;
}

uint32_t vector_from_stream(vector_t* vector_ptr, FILE* stream, char* name, uint32_t capacity, uint32_t size) {
    if (size > capacity) {
        fprintf(stderr, "Requested size exceeds capacity in vector_from_stream for %s\n", name);
        return 0;
    }

    if (!vector_init(vector_ptr, capacity)) {
        fprintf(stderr, "Failed to init vector %s in vector_from_stream\n", name);
        return 0;
    }

    if (stream == stdin) {
        printf("Enter %u values for vector %s: ", size, name);
        fflush(stdout);
    }

    for (uint32_t i = 0; i < size; i++) {
        if (fscanf(stream, "%lf", vector_ptr->data + i) != 1) {
            fprintf(stderr, "Failed to read element %u of vector %s\n", i, name);
            vector_free(vector_ptr);
            return 0;
        }
    }

    return 1;
}

uint32_t vector_as_gsl_view(vector_t* vector_ptr, uint32_t offset, uint32_t length, gsl_vector_view* view_ptr) {
    if (!vector_ptr || !view_ptr || !vector_ptr->data) {
        return 0;
    }

    if (offset + length > vector_ptr->capacity) {
        fprintf(stderr, "View exceeds capacity\n");
        return 0;
    }

    *view_ptr = gsl_vector_view_array(vector_ptr->data + offset, length);
    return 1;
}

double vector_get(const vector_t* vector_ptr, uint32_t i) {
    if (vector_ptr && i < vector_ptr->capacity && vector_ptr->data) {
        return vector_ptr->data[i];
    }

    __builtin_trap();
    return 0.0;
}

uint32_t vector_set(vector_t* vector_ptr, uint32_t i, double value) {
    if (!vector_ptr || i >= vector_ptr->capacity || !vector_ptr->data) {
        return 0;
    }

    vector_ptr->data[i] = value;
    return 1;
}

void vector_free(vector_t* vector_ptr) {
    if (!vector_ptr) {
        return;
    }

    free(vector_ptr->data);
    vector_ptr->data = NULL;
    vector_ptr->capacity = 0;
}