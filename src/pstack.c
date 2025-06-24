#include "pstack.h"
#include <stdio.h>
#include <stdlib.h>

/* PLIST */
struct pnode {
    problem_t p;
    plist_t next;
};

plist_t plist_new(const problem_t problem) {
    if (!problem) {
        return NULL;
    }

    plist_t l = (plist_t)malloc(sizeof(_pnode));
    if (!l) {
        return NULL;
    }

    l->p = problem_duplicate(problem);
    l->next = NULL;

    return l;
}

uint32_t plist_empty(const plist_t l) {
    return !l;
}

uint32_t plist_insert(plist_t* lp, plist_t n) {
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

plist_t plist_remove(plist_t* lp) {
    if (!lp || plist_empty(*lp)) {
        return NULL;
    }

    plist_t removed = *lp;
    *lp = (*lp)->next;
    removed->next = NULL;

    return removed;
}

void plist_free(plist_t* lp) {
    if (!lp || !*lp) {
        return;
    }

    while (*lp) {
        plist_t next = (*lp)->next;
        problem_free(&(*lp)->p);
        free(*lp);
        *lp = next;
    }
}

/* PSTACK */
struct pstack {
    plist_t top;
    uint32_t size;
};

pstack_t pstack_new() {
    pstack_t s = (pstack_t)malloc(sizeof(_pstack));
    if (!s) {
        return NULL;
    }

    s->top = NULL;
    s->size = 0;

    return s;
}

uint32_t pstack_empty(const pstack_t s) {
    return !s || plist_empty(s->top);
}

uint32_t pstack_push(const pstack_t s, const problem_t problem) {
    if (!s) {
        return 0;
    }

    plist_t n = plist_new(problem);
    if (!n) {
        return 0;
    }

    if (plist_insert(&s->top, n)) {
        s->size++;
        return 1;
    }

    free(n);
    return 0;
}

problem_t pstack_pop(const pstack_t s) {
    if (!s) {
        return NULL;
    }

    plist_t removed = plist_remove(&s->top);
    if (!removed) {
        return NULL;
    }

    problem_t p = removed->p;
    free(removed);

    s->size--;

    return p;
}

uint32_t pstack_size(const pstack_t s) {
    if (!s || pstack_empty(s)) {
        return 0;
    } else {
        return s->size;
    }
}

void pstack_free(pstack_t* sp) {
    if (!sp) {
        return;
    }

    plist_free(&(*sp)->top);
    free(*sp);
    *sp = NULL;
}