#ifndef BB_STACK_H
#define BB_STACK_H

#include <branch_bound/node.h>

/* PLIST */
typedef struct pnode {
    bb_node_t bb_node;
    struct pnode* next;
} pnode_t;

pnode_t* pnode_new(bb_node_t bb_node);
uint32_t plist_empty(const pnode_t* l);
uint32_t plist_insert(pnode_t** lp, pnode_t* n);
pnode_t* plist_remove(pnode_t** lp);
void plist_free(pnode_t** lp);

/* PSTACK */
typedef struct pstack {
    pnode_t* top;
    uint32_t size;
} pstack_t;

uint32_t pstack_init(pstack_t* pstack_ptr);
uint32_t pstack_empty(const pstack_t* pstack_ptr);
uint32_t pstack_push(pstack_t* pstack_ptr, bb_node_t bb_node);
uint32_t pstack_pop(pstack_t* pstack_ptr, bb_node_t* bb_node_ptr);
uint32_t pstack_size(const pstack_t* pstack_ptr);
void pstack_free(pstack_t* pstack_ptr);

#endif