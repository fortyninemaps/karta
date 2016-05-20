
cdef extern from "linkedlist.c":
    cdef struct LinkedListNode:
        LinkedListNode *next
        void *value
    LinkedListNode *ll_new(int size)
    LinkedListNode *ll_last(LinkedListNode* lln)
    void ll_append(LinkedListNode *lln, void *ptr)
    LinkedListNode *ll_destroy_head(LinkedListNode *lln)
    void ll_destroy(LinkedListNode *lln)

