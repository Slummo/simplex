#ifndef RC_H
#define RC_H

typedef void (*free_fn)(void* data);
typedef struct rc rc_t;

/// @brief Creates a new reference-counting pointer
/// @param data Raw pointer to the data
/// @param drop A function to free the data
/// @return A reference-counting pointer of data on success, NULL on fail
rc_t* rc_new(void* data, free_fn drop);

/// @brief Clones a reference-counting pointer
/// @param rc Refence-counting pointer
/// @return A cloned pointer or NULL if rc is NULL
rc_t* rc_clone(rc_t* rc);

/// @brief Gets the raw data contained in rc
/// @param rc Reference-counting pointer
/// @return The raw data on success, NULL if rc is NULL
const void* rc_data(const rc_t* rc);

/// @brief  Gets the raw data as mutable contained in rc
/// @param rc Reference-counting pointer
/// @return The raw mutable data on success, NULL if rc is NULL
void* rc_data_mut(rc_t* rc);

/// @brief Frees the reference-counting pointer, setting it to NULL
/// @param rc_ptr A raw pointer to the reference-counting pointer
void rc_free(rc_t** rc_ptr);

#endif