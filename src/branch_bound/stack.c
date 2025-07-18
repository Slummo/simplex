#include "branch_bound/stack.h"

#include <stdio.h>
#include <stdlib.h>

pnode_t* pnode_new(bb_node_t* bb_node_ptr) {
    if (!bb_node_ptr) {
        return NULL;
    }

    pnode_t* n = (pnode_t*)malloc(sizeof(pnode_t));
    if (!n) {
        return NULL;
    }

    n->bb_node_ptr = bb_node_ptr;
    n->next = NULL;

    return n;
}

uint32_t plist_empty(const pnode_t* l) {
    return !l;
}

uint32_t plist_insert(pnode_t** lp, pnode_t* n) {
    if (!lp) {
        return 0;
    }

    if (plist_empty(*lp)) {
        *lp = n;
    } else {
        n->next = *lp;
        *lp = n;
    }

    return 1;
}

pnode_t* plist_remove(pnode_t** lp) {
    if (!lp || plist_empty(*lp)) {
        return NULL;
    }

    pnode_t* removed = *lp;
    *lp = (*lp)->next;
    removed->next = NULL;

    return removed;
}

void plist_free(pnode_t** lp) {
    if (!lp || !*lp) {
        return;
    }

    while (*lp) {
        pnode_t* next = (*lp)->next;
        free(*lp);
        *lp = next;
    }
}

uint32_t pstack_init(pstack_t* pstack_ptr) {
    if (!pstack_ptr) {
        return 0;
    }

    pstack_ptr->top = NULL;
    pstack_ptr->size = 0;

    return 1;
}

uint32_t pstack_empty(const pstack_t* pstack_ptr) {
    return !pstack_ptr || plist_empty(pstack_ptr->top);
}

uint32_t pstack_push(pstack_t* pstack_ptr, bb_node_t* bb_node_ptr) {
    if (!pstack_ptr || !bb_node_ptr) {
        return 0;
    }

    pnode_t* n = pnode_new(bb_node_ptr);
    if (!n) {
        return 0;
    }

    if (plist_insert(&pstack_ptr->top, n)) {
        pstack_ptr->size++;
        return 1;
    }

    return 0;
}

bb_node_t* pstack_pop(pstack_t* pstack_ptr) {
    if (!pstack_ptr) {
        return NULL;
    }

    pnode_t* removed = plist_remove(&pstack_ptr->top);
    if (!removed) {
        return NULL;
    }

    bb_node_t* bb_node_ptr = removed->bb_node_ptr;
    free(removed);

    pstack_ptr->size--;

    return bb_node_ptr;
}

uint32_t pstack_size(const pstack_t* pstack_ptr) {
    if (pstack_empty(pstack_ptr)) {
        return 0;
    } else {
        return pstack_ptr->size;
    }
}

void pstack_free(pstack_t* pstack_ptr) {
    if (!pstack_ptr) {
        return;
    }

    plist_free(&pstack_ptr->top);
    pstack_ptr->size = 0;
}