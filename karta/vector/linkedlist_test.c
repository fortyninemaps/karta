#include <stdio.h>
#include "linkedlist.h"

int main() {
    char *string = "apple";
    printf("creating list with a string\n  %s\n", string);
    LLNode *list = ll_new(string);

    int number = 42;
    printf("appending an integer\n  %d\n", number);
    ll_append(list, &number);

    printf("printing value of first item\n");
    printf("  %s\n", (char*) list->value);

    printf("printing value of last item\n");
    printf("  %d\n", *((int*) ll_last(list)->value));

    printf("cleaning up\n");
    ll_destroy_head(list);
    return 0;
}
