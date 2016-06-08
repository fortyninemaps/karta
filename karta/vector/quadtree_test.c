#include <stdio.h>
#include <stdlib.h>
#include "quadtree.h"


int main() {

    srand(49);

    Bbox *bb = qt_new_bbox(0, 0, 1, 1);
    NodePtrUnion root;
    NonleafNode *retnode;
    root.leafnode = qt_new_leaf(50, bb);

    int i;
    Position p;
    for (i=0; i!=1000000; i++) {
        Position *p = (Position*) malloc(sizeof(Position));
        p->id = i;
        p->x = rand()/(double) RAND_MAX;
        p->y = rand()/(double) RAND_MAX;

        retnode = qt_insert(root, *p);
        if (retnode != NULL) {
            root.nonleafnode = retnode;
        }

        free(p);
    }

    qt_free_node(root);
    return 0;
}
