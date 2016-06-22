// pool data structure
#include <stdlib.h>
#include <stdio.h>

typedef struct PointerPool {
    int size;
    char** members;
    int count;
} Pool;


Pool* pool_new(int membersize, int size) {
    Pool *pool = malloc(sizeof(Pool));
    pool->size = size;
    pool->count = 0;
    pool->members = malloc(size*membersize);
    return pool;
}

int pool_resize(Pool *pool, int membersize) {
    char **newmembers = malloc((pool->size)*2 * membersize);
    int i = 0;
    while (i != pool->size) {
        newmembers[i] = pool->members[i];
        i += 1;
    }
    free(pool->members);
    pool->members = newmembers;
    pool->size *= 2;
    return 0;
}

int pool_add(Pool *pool, char *item) {
    int err = 0;
    if (pool->size == pool->count) {
        err = pool_resize(pool, sizeof(item));
    }
    if (err != 0) {
        return err;
    } else {
        pool->members[pool->count] = item;
        pool->count += 1;
        return 0;
    }
}

// remove item from pool and return reference
char* pool_pop(Pool *pool, int index) {
    if (index >= pool->count) {
        printf("pool error: pop index beyond pool size\n");
        exit(1);
    }
    char *ptr = pool->members[index];
    if (index != pool->count-1) {
        // swap final item into place
        pool->members[index] = pool->members[pool->count-1];
    }
    pool->count -= 1;
    return ptr;
}

// get reference without removing from pool
char* pool_get(Pool *pool, int index) {
    if (index >= pool->count) {
        printf("pool error: get index beyond pool size\n");
        exit(1);
    }
    char *ptr = pool->members[index];
    return ptr;
}

void pool_destroy(Pool *pool) {
    free(pool->members);
    free(pool);
}

