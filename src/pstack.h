#ifndef PSTACK_H
#define PSTACK_H

#include "problem.h"
#include <stdint.h>

/* PLIST */
typedef struct pnode _pnode, *plist_t;

plist_t plist_new(const problem_t problem);
uint32_t plist_empty(const plist_t l);
uint32_t plist_insert(plist_t* lp, plist_t n);
plist_t plist_remove(plist_t* lp);
void plist_free(plist_t* lp);

/* PSTACK */
typedef struct pstack _pstack, *pstack_t;

pstack_t pstack_new();
uint32_t pstack_empty(const pstack_t s);
uint32_t pstack_push(const pstack_t s, const problem_t problem);
problem_t pstack_pop(const pstack_t s);
uint32_t pstack_size(const pstack_t s);
void pstack_free(pstack_t* sp);

#endif