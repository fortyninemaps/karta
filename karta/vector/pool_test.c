#include <stdio.h>
#include "pool.h"

int main() {
    printf("starting\n");
    Pool *pool = pool_new(sizeof(double*), 2);
    printf("allocated\n");

    printf("count: %d\n", pool->count);
    printf("size:  %d\n", pool->size);

    int err = 0;
    double d1 = 1.0;
    err = pool_add(pool, (char*) &d1);
    printf("\terr: %d\n", err);
    double d2 = 2.0;
    err = pool_add(pool, (char*) &d2);
    printf("\terr: %d\n", err);
    double d3 = 3.0;
    err = pool_add(pool, (char*) &d3);
    printf("\terr: %d\n", err);

    printf("count: %d\n", pool->count);
    printf("size:  %d\n", pool->size);


    double r1, r2;
    r1 = *((double*) pool_pop(pool, 1));
    r2 = *((double*) pool_pop(pool, 1));
    printf("r1: %f\n", r1);
    printf("r2: %f\n", r2);

    printf("count: %d\n", pool->count);
    printf("size:  %d\n", pool->size);


    printf("clearing\n");
    pool_destroy(pool);
    return 0;
}

