#ifndef PSTACK_H
#define PSTACK_H

#include "problem.h"
#include <stdint.h>

/* PLIST */
typedef struct pnode pnode_t;

pnode_t* pnode_new(const problem_t* problem);
uint32_t plist_empty(const pnode_t* l);
uint32_t plist_insert(pnode_t** lp, pnode_t* n);
pnode_t* plist_remove(pnode_t** lp);
void plist_free(pnode_t** lp);

/* PSTACK */
typedef struct pstack pstack_t;

pstack_t* pstack_new();
uint32_t pstack_empty(const pstack_t* s);
uint32_t pstack_push(pstack_t* s, const problem_t* problem);
problem_t* pstack_pop(pstack_t* s);
uint32_t pstack_size(const pstack_t* s);
void pstack_free(pstack_t** sp);

#endif