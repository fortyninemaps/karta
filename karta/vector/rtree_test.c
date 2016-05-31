#include <stdio.h>
#include "rtree.h"

void print_indent(int indent) {
    for (int i = 0; i != indent; i++) {
        printf(" ");
    }
}

void print_tree(Node *tree, int indent) {
    if (tree->type == LEAF) {
        print_indent(indent);
        printf("leaf node with %d boxes\n", tree->count);
    } else {
        print_indent(indent);
        printf("nonleaf node with %d subnodes\n", tree->count);
        for (int i=0; i!=tree->count; i++) {
            print_tree((Node*) tree->children[i], indent+2);
        }
    }
}

int main() {
    float x = -10.0;
    float y;
    int i = 0;

    Node *root = rt_new_node(LEAF, LINEAR, 8, NULL);
    Bbox *bb = NULL;

    while (x < 10.0) {
        y = x*x - 2*x - 25.0;

        bb = rt_new_bbox();
        bb->xmin = x;
        bb->ymin = y;
        bb->xmax = x+0.1;
        bb->ymax = y+0.2;
        //printf("%p\n", bb);
        root = rt_insert(root, bb, i);

        x = x + 0.1;
        i++;
    }
    printf("added %d geometries\n", i);
    printf("summary:\n");
    print_tree(root, 2);
    rt_free(root);
    printf("memory free'd\n");
    return 0;
}
