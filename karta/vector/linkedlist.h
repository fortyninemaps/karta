// Simple linked list for use as a FIFO queue
// see crtree.pyx for use

#include <stdlib.h>

typedef struct LinkedListNode {
    struct LinkedListNode *next;
    void *value;
} LLNode;

// create a new linked list
static LLNode *ll_new(void *ptr) {
    LLNode *node;
    node = malloc(sizeof(LLNode));
    node->value = ptr;
    node->next = NULL;
    return node;
}

// get the last item of a linked list
static LLNode *ll_last(LLNode *lln) {
    LLNode *node;
    node = lln;
    while (node->next) {
        node = node->next;
    }
    return node;
}

// add an item to the end of a linked list
static void ll_append(LLNode *lln, void *ptr) {
    LLNode *child = ll_new(ptr);
    LLNode *last = ll_last(lln);
    last->next = child;
    child->next = NULL;
}

// free memory used by head of linked list and return second item or NULL
static LLNode *ll_destroy_head(LLNode *lln) {
    LLNode *next;
    if (lln->next) {
        next = lln->next;
    } else {
        next = NULL;
    }
    free(lln);
    return next;
}

// free memory used by LinkedListNode and every following LinkedListNode
static void ll_destroy(LLNode *lln) {
    if (lln->next) {
        ll_destroy(lln->next);
    }
    free(lln);
}

