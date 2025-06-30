#include "pstack.h"
#include <stdio.h>
#include <stdlib.h>

/* PLIST */
struct pnode {
    problem_t* p;
    pnode_t* next;
};

pnode_t* pnode_new(problem_t* problem) {
    if (!problem) {
        return NULL;
    }

    pnode_t* n = (pnode_t*)malloc(sizeof(pnode_t));
    if (!n) {
        return NULL;
    }

    n->p = problem;
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
        problem_free(&(*lp)->p);
        free(*lp);
        *lp = next;
    }
}

/* PSTACK */
struct pstack {
    pnode_t* top;
    uint32_t size;
};

pstack_t* pstack_new() {
    pstack_t* s = (pstack_t*)malloc(sizeof(pstack_t));
    if (!s) {
        return NULL;
    }

    s->top = NULL;
    s->size = 0;

    return s;
}

uint32_t pstack_empty(const pstack_t* s) {
    return !s || plist_empty(s->top);
}

uint32_t pstack_push(pstack_t* s, problem_t* problem) {
    if (!s) {
        return 0;
    }

    pnode_t* n = pnode_new(problem);
    if (!n) {
        return 0;
    }

    if (plist_insert(&s->top, n)) {
        s->size++;
        return 1;
    }

    return 0;
}

problem_t* pstack_pop(pstack_t* s) {
    if (!s) {
        return NULL;
    }

    pnode_t* removed = plist_remove(&s->top);
    if (!removed) {
        return NULL;
    }

    problem_t* p = removed->p;
    free(removed);

    s->size--;

    return p;
}

uint32_t pstack_size(const pstack_t* s) {
    if (!s || pstack_empty(s)) {
        return 0;
    } else {
        return s->size;
    }
}

void pstack_free(pstack_t** sp) {
    if (!sp) {
        return;
    }

    plist_free(&(*sp)->top);
    free(*sp);
    *sp = NULL;
}