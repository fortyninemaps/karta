#include <stdio.h>
#include <stdlib.h>
#include "quadtree.h"


int main() {

    srand(49);

    Bbox *bb = qt_new_bbox(0, 0, 1, 1);
    NodePtrUnion root;
    NonleafNode *retnode;
    root.leafnode = qt_new_leaf(5, bb);
    NodeType type = LEAF;

    //printf("sizeof leaf: %u\n", sizeof *root.leafnode);

    int i;
    double x, y;
    Position p;
    for (i=0; i!=50; i++) {
        Position *p = (Position*) malloc(sizeof(Position));
        p->x = rand()/(double) RAND_MAX;
        p->y = rand()/(double) RAND_MAX;
        //printf("%.3f %.3f %p\n", p->x, p->y, p);

        retnode = qt_insert(root, *p);
        if ((type == LEAF) && (retnode != NULL)) {
            printf("split root (should only happen once)\n");
            type = NONLEAF;
            root.nonleafnode = retnode;
        }

        free(p);
    }

    //printf("sizeof nonleaf: %u\n", sizeof *root.nonleafnode);

    qt_free_node(root);
    return 0;
}
