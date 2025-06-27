#ifndef RC_H
#define RC_H

typedef void (*free_fn)(void* data);
typedef struct rc _rc, *rc_t;

rc_t rc_new(void* data, free_fn drop);
rc_t rc_clone(rc_t rc);
const void* rc_data(rc_t rc);
void* rc_data_mut(rc_t rc);
void rc_free(rc_t* rc_ptr);

#endif