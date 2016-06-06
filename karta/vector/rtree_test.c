#include <stdio.h>
#include "rtree.h"

void print_indent(int indent) {
    for (int i = 0; i != indent; i++) {
        printf(" ");
    }
}

void print_tree(Node *tree, int indent) {
    int i;
    if (tree->type == LEAF) {
        print_indent(indent);
        printf("leaf node with %d boxes (%p, parent is %p)\n", tree->count, tree, tree->parent);
        for (i=0; i!=tree->count; i++) {
            print_indent(indent+2);
            printf("%d %p\n", tree->indices[i], ((Bbox**) tree->children)[i]);
            // print_indent(indent+4);
            // print_bbox(((Bbox**) tree->children)[i]);
        }
    } else {
        print_indent(indent);
        printf("nonleaf node with %d subnodes (%p, parent is %p)\n", tree->count, tree, tree->parent);
        for (i=0; i!=tree->count; i++) {
            print_tree((Node*) tree->children[i], indent+2);
        }
    }
}

int test_is_within() {
    Bbox *bb = rt_new_bbox();
    Bbox *bb1 = rt_new_bbox();
    Bbox *bb2 = rt_new_bbox();

    bb->xmin = -10;
    bb->xmax = 10;
    bb->ymin = -10;
    bb->ymax = 10;

    bb1->xmin = -5;
    bb1->xmax = 5;
    bb1->ymin = -5;
    bb1->ymax = 5;

    bb2->xmin = -11;
    bb2->xmax = 5;
    bb2->ymin = -5;
    bb2->ymax = 5;

    if (!is_within(bb, bb1)) {
        printf("fail is_within test 1\n");
        return 1;
    }
    if (is_within(bb, bb2)) {
        printf("fail is_within test 2\n");
        return 1;
    }
    free(bb);
    free(bb1);
    free(bb2);
    return 0;
}

int test_is_overlapping() {
    Bbox *bb = rt_new_bbox();
    Bbox *bb1 = rt_new_bbox();
    Bbox *bb2 = rt_new_bbox();
    Bbox *bb3 = rt_new_bbox();

    bb->xmin = -10;
    bb->xmax = 10;
    bb->ymin = -10;
    bb->ymax = 10;

    bb1->xmin = -5;
    bb1->xmax = 5;
    bb1->ymin = -5;
    bb1->ymax = 5;

    bb2->xmin = 0;
    bb2->xmax = 15;
    bb2->ymin = 0;
    bb2->ymax = 15;

    bb3->xmin = -20;
    bb3->xmax = -15;
    bb3->ymin = 8;
    bb3->ymax = 12;

    if (!is_overlapping(bb, bb1)) {
        printf("fail is_overlapping test 1\n");
        return 1;
    }
    if (!is_overlapping(bb, bb2)) {
        printf("fail is_overlapping test 2\n");
        return 1;
    }
    if (is_overlapping(bb, bb3)) {
        printf("fail is_overlapping test 3\n");
        return 1;
    }
    free(bb);
    free(bb1);
    free(bb2);
    free(bb3);
    return 0;
}

int count_geometries(Node *node) {
    Pool *pool = pool_new(sizeof(Node*), 4);
    Node *active_node = NULL;
    int count = 0;
    pool_add(pool, (char*) node);
    while (pool->count != 0) {

        active_node = (Node*) pool_pop(pool, pool->count-1);
        if (active_node->type == LEAF) {
            count = count + active_node->count;
        } else if (active_node->type == NONLEAF) {
            for (int i=0; i!=active_node->count; i++) {
                pool_add(pool, (char*) ((Node**) active_node->children)[i]);
            }
        }
    }
    pool_destroy(pool);
    return count;
}

int test_tree_construction() {
    float x = -10.0;
    float y;
    int i = 0;

    Node *root = rt_new_node(LEAF, LINEAR, 25, NULL);
    Bbox *bb = NULL;

    while (x < 10.0) {
        y = x*x;

        bb = rt_new_bbox();
        bb->xmin = x;
        bb->ymin = y;
        bb->xmax = x+0.1;
        bb->ymax = y+0.1;
        root = rt_insert(root, bb, i);

        x = x + 0.2;
        i++;
    }
    if (i != count_geometries(root)) {
        printf("failed geometry count\n");
        rt_free(root);
        return 1;
    } else {
        rt_free(root);
        return 0;
    }
}

int main() {
    int ip;
    ip = test_is_within();
    if (ip != 0) {
        printf("FAILED IS_WITHIN TEST\n");
        exit(1);
    }

    ip = test_is_overlapping();
    if (ip != 0) {
        printf("FAILED IS_OVERLAPPING TEST\n");
        exit(1);
    }

    ip = test_tree_construction();
    if (ip != 0) {
        printf("FAILED TREE_CONSTRUCTION TEST\n");
        exit(1);
    }

    // Tree search test
    float x = -10.0;
    float y;
    int i = 0;

    Node *root = rt_new_node(LEAF, LINEAR, 17, NULL);
    Bbox *bb = NULL;

    while (x < 10.0) {
        y = x*x;

        bb = rt_new_bbox();
        bb->xmin = x;
        bb->ymin = y;
        bb->xmax = x+0.1;
        bb->ymax = y+0.1;
        //printf("adding %p\n", bb);
        root = rt_insert(root, bb, i);

        x = x + 0.1;
        i++;
    }
    printf("added %d geometries\n", i);
    printf("(tree contains %d geometries)\n", count_geometries(root));

    //printf("summary:\n");
    //print_tree(root, 2);
    print_bbox(root->bbox);

    printf("searching\n");
    bb = rt_new_bbox();
    bb->xmin = -10;
    bb->xmax = 0;
    bb->ymin = 0;
    bb->ymax = 30;
    Pool *results = rt_search_within(root, bb, -1);
    printf("number of results: %d\n", results->count);

    free(bb);
    pool_destroy(results);
    rt_free(root);
    printf("memory free'd\n");
    return 0;
}
