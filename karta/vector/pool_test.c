#include <stdlib.h>
#include <stdio.h>
#include "pool.h"

int main() {
    printf("starting\n");
    Pool *pool = pool_new(sizeof(double*), 2);
    printf("allocated\n");

    printf("count: %d\n", pool->count);
    printf("size:  %d\n", pool->size);

    int i, err;
    double a[10];
    for (i=0; i!=10; i++) {
        a[i] = i*1.0;
        printf("%f\n", a[i]);
        err = pool_add(pool, (char*) &a[i]);
        if (err!=0) {
            printf("\tALLOCATION ERROR: %d\n", err);
            exit(1);
        }
    }

    printf("count: %d\n", pool->count);
    printf("size:  %d\n", pool->size);

    double r1, r2, r3;
    for (i=0; i!=pool->count; i++) {
        printf("item %d: %f\n", i, *((double**) pool->members)[i]);
    }
    r1 = *((double*) pool_pop(pool, 1));
    printf("r1: %f\n", r1);
    for (i=0; i!=pool->count; i++) {
        printf("item %d: %f\n", i, *((double**) pool->members)[i]);
    }
    r2 = *((double*) pool_pop(pool, 2));
    printf("r2: %f\n", r2);
    for (i=0; i!=pool->count; i++) {
        printf("item %d: %f\n", i, *((double**) pool->members)[i]);
    }

    r3 = *((double*) pool_pop(pool, 1));
    printf("r3: %f\n", r3);
    for (i=0; i!=pool->count; i++) {
        printf("item %d: %f\n", i, *((double**) pool->members)[i]);
    }

    printf("count: %d\n", pool->count);
    printf("size:  %d\n", pool->size);

    printf("clearing pool\n");
    pool_destroy(pool);
    return 0;
}

