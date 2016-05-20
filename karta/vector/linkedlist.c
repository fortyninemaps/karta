// Simple linked list for use as a FIFO queue
// see crtree.pyx for use

#include <stdlib.h>

struct LinkedListNode {
    struct LinkedListNode* next;
    void *value;
};

static struct LinkedListNode *ll_new(int size) {
    struct LinkedListNode *node;
    node = (struct LinkedListNode*) malloc(sizeof(struct LinkedListNode));
    node->next = NULL;
    return node;
}

static struct LinkedListNode *ll_last(struct LinkedListNode* lln) {
    struct LinkedListNode *node;
    node = lln;
    while (node->next) {
        node = node->next;
    }
    return node;
}

static void ll_append(struct LinkedListNode* lln, void* ptr) {
    struct LinkedListNode *child = ll_new(sizeof(ptr));
    struct LinkedListNode *last = ll_last(lln);
    child->value = ptr;
    child->next = NULL;
    last->next = child;
}

static struct LinkedListNode *ll_destroy_head(struct LinkedListNode* lln) {
    struct LinkedListNode *next;
    if (lln->next) {
        next = lln->next;
    } else {
        next = NULL;
    }
    free(lln);
    return next;
}

static void ll_destroy(struct LinkedListNode* lln) {
    if (lln->next) {
        ll_destroy(lln->next);
    }
    free(lln);
}

