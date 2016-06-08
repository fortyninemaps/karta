// R-tree implementation for indexing Karta geometries

#include <stdlib.h>
#include <stdio.h>
#include "linkedlist.h"
#include "pool.h"

typedef enum { LINEAR, QUADRATIC } Strategy;

typedef enum { LEAF, NONLEAF } NodeType;

typedef struct BoundingBox {
    float xmin;
    float ymin;
    float xmax;
    float ymax;
} Bbox;

typedef struct RTreeNode {
    NodeType type;
    Strategy strategy;
    int count;
    int maxchildren;
    struct RTreeNode* parent;
    char** children;
    int* indices;
    Bbox* bbox;
} Node;

union ChildUnionTag {
    Node* node;
    Bbox* bbox;
} ChildUnion;

float minf(float a, float b) {
    return (a >= b) ? b : a;
}

float maxf(float a, float b) {
    return (a >= b) ? a : b;
}

// forward declarations
Node* rt_insert(Node*, Bbox*, int);
Node* rt_choose_leaf(Node*, Bbox*);
Node* rt_split(Node*);
Node* rt_split_nonleaf(Node*);
Node* rt_split_leaf(Node*);
int rt_adjust_tree(Node*, Node*, Node**, Node**);
float volume_expanded(Bbox*, Bbox*);
int linear_pick_seeds(int, Bbox**, int*, int*);
int linear_pick_next(int, Bbox**, Node*, Node*, int*, int*);
int is_within(Bbox*, Bbox*);
int is_overlapping(Bbox*, Bbox*);

// new_node allocates a new node struct
Node* rt_new_node(NodeType type, Strategy strategy, int maxchildren, Node* parent) {
    if (maxchildren<4) {
        printf("rtree error: maxchildren must be greater than 4");
        exit(1);
    }
    Node *node = malloc(sizeof(Node));
    node->type = type;
    node->strategy = strategy;
    node->count = 0;
    node->maxchildren = maxchildren;
    node->parent = parent;
    node->children = malloc((maxchildren+1) * sizeof(ChildUnion));
    node->indices = malloc((maxchildren+1) * sizeof(int));

    // TODO: use NaN values
    node->bbox = malloc(sizeof(Bbox));
    node->bbox->xmin = 1e38;
    node->bbox->xmax = -1e38;
    node->bbox->ymin = 1e38;
    node->bbox->ymax = -1e38;
    return node;
}

Bbox* rt_new_bbox() {
    Bbox* bb = malloc(sizeof(Bbox));
    return bb;
}

void rt_free(Node *node) {
    int i;
    if (node->type == LEAF) {
        for (i=0; i!=node->count; i++) {
            free(((Bbox**) node->children)[i]);
        }
        free(node->indices);
        free(node->bbox);
        free(node);
    } else if (node->type == NONLEAF) {
        for (i=0; i!=node->count; i++) {
            rt_free(((Node**) node->children)[i]);
        }
    }
}

void print_bbox(Bbox *bb) {
    printf("%.3f %.3f %.3f %.3f\n", bb->xmin, bb->ymin, bb->xmax, bb->ymax);
}

// destroy_node frees the memory associated with a single node
void rt_free_single(Node *node) {
    free(node->children);
    free(node->indices);
    free(node);
}

// add_node appends a node to the children of parent
// the return value is 1 if the resulting parent is overfilled
int add_node(Node *parent, Node *child) {
    int ret;
    parent->children[parent->count] = (char*) child;
    parent->count++;
    child->parent = parent;
    if (parent->count > parent->maxchildren) {
        ret = 1;
    } else {
        ret = 0;
    }
    return ret;
}

// search for geometries within a bbox
Pool *rt_search_within(Node *node, Bbox *bbox, int max_results) {
    Pool *index_pool = pool_new(sizeof(int*), 16);
    Pool *node_pool = pool_new(sizeof(Node*), 16);
    pool_add(node_pool, (char*) node);
    Node *active_node = NULL;
    Node *child_node = NULL;
    Bbox *child_bbox = NULL;
    int nfound = 0;
    int j;
    while (nfound != max_results) {
        active_node = (Node*) pool_pop(node_pool, node_pool->count-1);
        if (active_node->type == LEAF) {
            for (j=0; j!=active_node->count; j++){
                child_bbox = ((Bbox**) active_node->children)[j];
                if (is_within(bbox, child_bbox) == 1) {
                    pool_add(index_pool, (char*) &(active_node->indices[j]));
                    nfound++;
                }
            }
        } else if (active_node->type == NONLEAF) {
            for (j=0; j!=active_node->count; j++) {
                child_node = ((Node**) active_node->children)[j];
                if (is_overlapping(bbox, child_node->bbox) == 1) {
                    pool_add(node_pool, (char*) child_node);
                }
            }
        }

        if (node_pool->count == 0) {
            break;
        }
    }
    pool_destroy(node_pool);
    return index_pool;
}

// search for geometries overlapping a bbox
Pool *rt_search_overlapping(Node *node, Bbox *bbox, int max_results) {
    Pool *index_pool = pool_new(sizeof(int*), 16);
    Pool *node_pool = pool_new(sizeof(Node*), 16);
    pool_add(node_pool, (char*) node);
    Node *active_node = NULL;
    Node *child_node = NULL;
    Bbox *child_bbox = NULL;
    int i = 0;
    int j;
    while (i != max_results) {
        active_node = (Node*) pool_pop(node_pool, node_pool->count-1);
        if (active_node->type == LEAF) {
            for (j=0; j!=active_node->count; j++){
                child_bbox = ((Bbox**) active_node->children)[j];
                if (is_overlapping(bbox, child_bbox) == 1) {
                    pool_add(index_pool, (char*) &(active_node->indices[j]));
                }
            }
        } else if (active_node->type == NONLEAF) {
            for (j=0; j!=active_node->count; j++) {
                child_node = ((Node**) active_node->children)[j];
                if (is_overlapping(bbox, child_node->bbox) == 1) {
                    pool_add(node_pool, (char*) child_node);
                }
            }
        }

        if (node_pool->count == 0) {
            break;
        }
    }
    pool_destroy(node_pool);
    return index_pool;
}

// add a bbox to rtree, returning a reference to the root node
Node* rt_insert(Node *root, Bbox *bbox, int index) {
    Node *sibling = NULL;
    Node *leaf = rt_choose_leaf(root, bbox);
    Node **outnode = malloc(sizeof(Node*));
    Node **outsibling = malloc(sizeof(Node*));
    Node *returned_root;

    leaf->children[leaf->count] = (char*) bbox;
    leaf->indices[leaf->count] = index;
    leaf->count++;

    leaf->bbox->xmin = minf(leaf->bbox->xmin, bbox->xmin);
    leaf->bbox->xmax = maxf(leaf->bbox->xmax, bbox->xmax);
    leaf->bbox->ymin = minf(leaf->bbox->ymin, bbox->ymin);
    leaf->bbox->ymax = maxf(leaf->bbox->ymax, bbox->ymax);

    if (leaf->count > leaf->maxchildren) {
        sibling = rt_split(leaf);
    }

    int ret = rt_adjust_tree(leaf, sibling, outnode, outsibling);
    if (ret == 1) {
        returned_root = rt_new_node(NONLEAF, root->strategy, root->maxchildren, NULL);
        add_node(returned_root, *outnode);
        add_node(returned_root, *outsibling);
    } else {
        returned_root = *outnode;
    }
    free(outnode);
    free(outsibling);
    return returned_root;
}

// chose a leaf node below node in which to place bbox
Node* rt_choose_leaf(Node *node, Bbox *bbox) {
    int i, iminexpanded;
    float vol, cur_vol;
    while (node->type != LEAF) {
        iminexpanded = 0;
        vol = volume_expanded(((Node**) node->children)[0]->bbox, bbox);
        for (i=1; i!=node->count; i++) {
            cur_vol = volume_expanded(((Node**) node->children)[i]->bbox, bbox);
            if (cur_vol < vol) {
                iminexpanded = i;
                vol = cur_vol;
            }
        }
        node = ((Node**) node->children)[iminexpanded];
    }
    return node;
}

// propagate changes upward. return value is 1 if the root node was split
//
// after returning, **outnode is the tree root, and **outsibling is NULL or a
// root sibling if the root was split
int rt_adjust_tree(Node *node, Node *sibling, Node **outnode, Node **outsibling) {
    int i, ret = 0;
    int addretval = 0;
    float xmin, xmax, ymin, ymax;
    Node *parent;
    while (node->parent != NULL) { // while not the root node
        // make sure the parent's bbox is up to date
        parent = node->parent;
        if (sibling == NULL) {
            xmin = ((Node**) parent->children)[0]->bbox->xmin;
            xmax = ((Node**) parent->children)[0]->bbox->xmax;
            ymin = ((Node**) parent->children)[0]->bbox->ymin;
            ymax = ((Node**) parent->children)[0]->bbox->ymax;
            i = 1;
        } else {
            xmin = sibling->bbox->xmin;
            xmax = sibling->bbox->xmax;
            ymin = sibling->bbox->ymin;
            ymax = sibling->bbox->ymax;
            i = 0;
        }
        while (i != parent->count) {
            xmin = minf(xmin, ((Node**) parent->children)[i]->bbox->xmin);
            xmax = maxf(xmax, ((Node**) parent->children)[i]->bbox->xmax);
            ymin = minf(ymin, ((Node**) parent->children)[i]->bbox->ymin);
            ymax = maxf(ymax, ((Node**) parent->children)[i]->bbox->ymax);
            i++;
        }
        parent->bbox->xmin = xmin;
        parent->bbox->xmax = xmax;
        parent->bbox->ymin = ymin;
        parent->bbox->ymax = ymax;

        // if sibling provided, add it to parent
        if (sibling != NULL) {
            addretval = add_node(node->parent, sibling);
            if (addretval == 1) {
                sibling = rt_split(node->parent);
            } else {
                sibling = NULL;
            }
        }

        // propagate upward
        node = node->parent;
    }

    if (sibling != NULL) {
        ret = 1;
    }

    *outnode = node;
    *outsibling = sibling;
    return ret;
}

int rt_tighten_bbox(Node *node) {
    if (node->count == 0) {
        return 1;
    }
    int i;
    if (node->type == LEAF) {
        node->bbox->xmin = ((Bbox**) node->children)[0]->xmin;
        node->bbox->xmax = ((Bbox**) node->children)[0]->xmax;
        node->bbox->ymin = ((Bbox**) node->children)[0]->ymin;
        node->bbox->ymax = ((Bbox**) node->children)[0]->ymax;
        for (i=1; i!=node->count; i++) {
            node->bbox->xmin = minf(node->bbox->xmin,
                                    ((Bbox*) node->children[i])->xmin);
            node->bbox->xmax = maxf(node->bbox->xmax,
                                    ((Bbox*) node->children[i])->xmax);
            node->bbox->ymin = minf(node->bbox->ymin,
                                    ((Bbox*) node->children[i])->ymin);
            node->bbox->ymax = maxf(node->bbox->ymax,
                                    ((Bbox*) node->children[i])->ymax);
        }
    } else if (node->type == NONLEAF) {
        node->bbox->xmin = ((Node**) node->children)[0]->bbox->xmin;
        node->bbox->xmax = ((Node**) node->children)[0]->bbox->xmax;
        node->bbox->ymin = ((Node**) node->children)[0]->bbox->ymin;
        node->bbox->ymax = ((Node**) node->children)[0]->bbox->ymax;
        for (i=1; i!=node->count; i++) {
            node->bbox->xmin = minf(node->bbox->xmin,
                                    ((Node**) node->children)[i]->bbox->xmin);
            node->bbox->xmax = maxf(node->bbox->xmax,
                                    ((Node**) node->children)[i]->bbox->xmax);
            node->bbox->ymin = minf(node->bbox->ymin,
                                    ((Node**) node->children)[i]->bbox->ymin);
            node->bbox->ymax = maxf(node->bbox->ymax,
                                    ((Node**) node->children)[i]->bbox->ymax);
        }

    }
    return 0;
}

// split an Rtree node, returning the new node. The original node and the new
// node combined contain the contents originally in node
Node* rt_split(Node *node) {
    Node *sibling = NULL;
    if (node->type == LEAF) {
        sibling = rt_split_leaf(node);
    } else if (node->type == NONLEAF) {
        sibling = rt_split_nonleaf(node);
    }
    return sibling;
}

// split_leaf redistributes children between a node and a new sibling node. The
// identity of the original node is not changed, although some of its children
// move to the sibling.
Node* rt_split_leaf(Node *node) {
    if (node->type != LEAF) {
        printf("rtree error: node passed to rt_split_leaf must be leaf\n");
        exit(1);
    }
    int i;
    // create a new leaf node to take some of node's children
    Node *sibling = rt_new_node(node->type, node->strategy, node->maxchildren, node->parent);

    // copy an array of bbox indices to back the pointer pool
    int *index_array;
    index_array = malloc((node->maxchildren+1) * sizeof(int));
    for (i=0; i!=node->maxchildren+1; i++) {
        index_array[i] = node->indices[i];
    }

    // construct a pool of the undistributed bounding boxes
    Pool *pool = pool_new(sizeof(Bbox*), node->maxchildren+1);
    Pool *idx_pool = pool_new(sizeof(int*), node->maxchildren+1);
    for (i=0; i!=pool->size; i++) {
        pool_add(pool, node->children[i]);
        pool_add(idx_pool, (char*) &index_array[i]);
    }
    node->count = 0;    // effectively empties node

    // select seeds for split
    int seed0 = -1, seed1 = -1;
    int err;
    err = linear_pick_seeds(pool->count, (Bbox**) pool->members, &seed0, &seed1);
    if (err != 0) {
        printf("warning: node geometries have degenerate dimension\n");
        seed0 = 0;
        seed1 = 1;
    }

    rt_insert(node, (Bbox*) pool->members[seed0], *((int*) idx_pool->members[seed0]));
    rt_insert(sibling, (Bbox*) pool->members[seed1], *((int*) idx_pool->members[seed1]));
    node->bbox->xmin = ((Bbox**) node->children)[0]->xmin;
    node->bbox->ymin = ((Bbox**) node->children)[0]->ymin;
    node->bbox->xmax = ((Bbox**) node->children)[0]->xmax;
    node->bbox->ymax = ((Bbox**) node->children)[0]->ymax;
    sibling->bbox->xmin = ((Bbox**) sibling->children)[0]->xmin;
    sibling->bbox->ymin = ((Bbox**) sibling->children)[0]->ymin;
    sibling->bbox->xmax = ((Bbox**) sibling->children)[0]->xmax;
    sibling->bbox->ymax = ((Bbox**) sibling->children)[0]->ymax;

    // remove seeds from pool in reverse order
    if (seed1 > seed0) {
        pool_pop(pool, seed1);
        pool_pop(pool, seed0);
        pool_pop(idx_pool, seed1);
        pool_pop(idx_pool, seed0);
    } else {
        pool_pop(pool, seed0);
        pool_pop(pool, seed1);
        pool_pop(idx_pool, seed0);
        pool_pop(idx_pool, seed1);
    }

    // choose bounding boxes from pool and place each in best of node, sibling
    int inextchild = -1;
    int ibestnode = -1;
    int *index;
    Node *target = NULL;
    Bbox *bbox = NULL;
    while (pool->count != 0) {
        err = linear_pick_next(pool->count, (Bbox**) pool->members, node, sibling,
                               &inextchild, &ibestnode);
        if (err != 0) {
            printf("rtree error: linear_pick_next\n");
            exit(err);
        }

        // copy the indicated bounding box from pool to the best node
        if (ibestnode == 0) {
            target = node;
        } else if (ibestnode == 1) {
            target = sibling;
        }

        bbox = (Bbox*) pool_pop(pool, inextchild);
        index = (int*) pool_pop(idx_pool, inextchild);
        target->children[target->count] = (char*) bbox;
        target->indices[target->count] = *index;
        target->count++;

        // after the index pool is popped, it's contents are correct
        // somehow at the next iteration, the final value is getting messed up

        target->bbox->xmin = minf(target->bbox->xmin, bbox->xmin);
        target->bbox->xmax = maxf(target->bbox->xmax, bbox->xmax);
        target->bbox->ymin = minf(target->bbox->ymin, bbox->ymin);
        target->bbox->ymax = maxf(target->bbox->ymax, bbox->ymax);

        // prevent either node from having fewer than 2 children
        if ((pool->count + node->count) == 2) {
            while (pool->count != 0) {
                node->children[node->count] = (char*) pool_pop(pool, pool->count-1);
                node->indices[node->count] = *((int*) pool_pop(idx_pool, idx_pool->count-1));
                node->count++;
            }
        } else if ((pool->count + sibling->count) == 2) {
            while (pool->count != 0) {
                sibling->children[sibling->count] = (char*) pool_pop(pool, pool->count-1);
                sibling->indices[sibling->count] = *((int*) pool_pop(idx_pool, idx_pool->count-1));
                sibling->count++;
            }
        }
    }

    free(index_array);
    pool_destroy(pool);
    pool_destroy(idx_pool);
    return sibling;
}

Node* rt_split_nonleaf(Node *node) {
    if (node->type != NONLEAF) {
        printf("rtree error: node passed to rt_split_nonleaf must be nonleaf\n");
        exit(1);
    }
    // create a new leaf node to take some of node's children
    Node *sibling = rt_new_node(node->type, node->strategy, node->maxchildren, node->parent);

    // construct a pool of the undistributed bounding boxes
    Pool *pool = pool_new(sizeof(Bbox*), node->maxchildren+1);
    Pool *node_pool = pool_new(sizeof(Node*), node->maxchildren+1);
    int i = 0;
    while (i != pool->size) {
        pool_add(pool, (char*) ((Node**) node->children)[i]->bbox);
        pool_add(node_pool, node->children[i]);
        i++;
    }
    node->count = 0;    // effectively empties node

    // select seeds for split
    int seed0 = -1, seed1 = -1;
    int err;
    err = linear_pick_seeds(pool->count, (Bbox**) pool->members, &seed0, &seed1);
    if (err != 0) {
        printf("warning: node geometries have degenerate dimension\n");
        seed0 = 0;
        seed1 = 1;
    }

    // copy seed children to node and sibling
    Node *tmpnode;

    tmpnode = ((Node**) node_pool->members)[seed0];
    node->children[0] = (char*) tmpnode;
    tmpnode->parent = node;
    node->count++;

    tmpnode = ((Node**) node_pool->members)[seed1];
    sibling->children[0] = (char*) tmpnode;
    tmpnode->parent = sibling;
    sibling->count++;

    node->bbox->xmin = ((Node**) node->children)[0]->bbox->xmin;
    node->bbox->ymin = ((Node**) node->children)[0]->bbox->ymin;
    node->bbox->xmax = ((Node**) node->children)[0]->bbox->xmax;
    node->bbox->ymax = ((Node**) node->children)[0]->bbox->ymax;
    sibling->bbox->xmin = ((Node**) sibling->children)[0]->bbox->xmin;
    sibling->bbox->ymin = ((Node**) sibling->children)[0]->bbox->ymin;
    sibling->bbox->xmax = ((Node**) sibling->children)[0]->bbox->xmax;
    sibling->bbox->ymax = ((Node**) sibling->children)[0]->bbox->ymax;

    // remove seeds from pool in reverse order
    if (seed1 > seed0) {
        pool_pop(pool, seed1);
        pool_pop(pool, seed0);
        pool_pop(node_pool, seed1);
        pool_pop(node_pool, seed0);
    } else {
        pool_pop(pool, seed0);
        pool_pop(pool, seed1);
        pool_pop(node_pool, seed0);
        pool_pop(node_pool, seed1);
    }

    // choose bounding boxes from pool and place each in best of node, sibling
    int inextchild = -1;
    int ibestnode = -1;
    Node *fosternode = NULL;
    while (node_pool->count != 0) {
        err = linear_pick_next(pool->count, (Bbox**) pool->members, node, sibling,
                               &inextchild, &ibestnode);
        if (err != 0) {
            printf("rtree error: linear_pick_next\n");
            exit(err);
        }

        // copy the indicated node from pool to the best node
        if (ibestnode == 0) {
            tmpnode = node;
        } else {
            tmpnode = sibling;
        }
        fosternode = (Node*) pool_pop(node_pool, inextchild);
        pool_pop(pool, inextchild);

        tmpnode->children[tmpnode->count] = (char*) fosternode;
        tmpnode->count++;
        tmpnode->bbox->xmin = minf(tmpnode->bbox->xmin, fosternode->bbox->xmin);
        tmpnode->bbox->xmax = maxf(tmpnode->bbox->xmax, fosternode->bbox->xmax);
        tmpnode->bbox->ymin = minf(tmpnode->bbox->ymin, fosternode->bbox->ymin);
        tmpnode->bbox->ymax = maxf(tmpnode->bbox->ymax, fosternode->bbox->ymax);

        // prevent either node from having fewer than 2 children
        if ((node_pool->count + node->count) == 2) {
            while (node_pool->count != 0) {
                node->children[node->count] = (char*) pool_pop(node_pool, node_pool->count-1);
                node->count++;
            }
        } else if ((node_pool->count + sibling->count) == 2) {
            while (node_pool->count != 0) {
                sibling->children[sibling->count] = (char*) pool_pop(node_pool, node_pool->count-1);
                sibling->count++;
            }
        }
    }

    pool_destroy(pool);
    pool_destroy(node_pool);
    return sibling;
}

Bbox* union_bbox(Bbox **bboxes, int nbboxes) {
    float xmin = bboxes[0]->xmin;
    float xmax = bboxes[0]->xmax;
    float ymin = bboxes[0]->ymin;
    float ymax = bboxes[0]->ymax;
    int i;
    for (i=0; i!=nbboxes; i++) {
        xmin = minf(xmin, bboxes[i]->xmin);
        xmax = maxf(xmax, bboxes[i]->xmax);
        ymin = minf(ymin, bboxes[i]->ymin);
        ymax = maxf(ymax, bboxes[i]->ymax);
    }
    Bbox* bb = rt_new_bbox();
    bb->xmin = xmin;
    bb->ymin = xmax;
    bb->xmax = ymin;
    bb->ymax = ymax;
    return bb;
}

// linear_pick_seeds selects two bboxes that should be used as seed bboxes at a
// new rtree node
int linear_pick_seeds(int nbboxes, Bbox **bboxes, int *seed0, int *seed1) {
    int iminxmax = 0, imaxxmin = 0;
    int iminymax = 0, imaxymin = 0;
    float minxmin = bboxes[0]->xmin;
    float maxxmax = bboxes[0]->xmax;
    float minymin = bboxes[0]->ymin;
    float maxymax = bboxes[0]->ymax;
    int i;
    for (i=1; i!=nbboxes; i++) {
        if (bboxes[i]->xmin > bboxes[imaxxmin]->xmin) {
            imaxxmin = i;
        } else if (bboxes[i]->xmax < bboxes[iminxmax]->xmax) {
            iminxmax = i;
        }
        if (bboxes[i]->ymin > bboxes[imaxymin]->ymin) {
            imaxymin = i;
        } else if (bboxes[i]->ymax < bboxes[iminymax]->ymax) {
            iminymax = i;
        }

        minxmin = minf(minxmin, bboxes[i]->xmin);
        maxxmax = maxf(maxxmax, bboxes[i]->xmax);
        minymin = minf(minymin, bboxes[i]->ymin);
        maxymax = maxf(maxymax, bboxes[i]->ymax);
    }
    if ((maxxmax == minxmin) || (maxymax == minymin)) {
        return 1;
    }
    float xrat = (bboxes[imaxxmin]->xmin - bboxes[iminxmax]->xmax) / (maxxmax - minxmin);
    float yrat = (bboxes[imaxymin]->ymin - bboxes[iminymax]->ymax) / (maxymax - minymin);
    if (xrat > yrat) {
        *seed0 = imaxxmin;
        *seed1 = iminxmax;
    } else {
        *seed0 = imaxymin;
        *seed1 = iminymax;
    }
    return 0;
}

// linear_pick_next chooses the next bbox from a list to allocate to a node, and
// chooses the node to allocate to
int linear_pick_next(int nbboxes, Bbox **bboxes, Node *node0, Node *node1,
                     int *inextchild, int *ibestnode) {
    *inextchild = nbboxes-1;
    Bbox *bbox = bboxes[nbboxes-1];
    float v0 = volume_expanded(node0->bbox, bbox);
    float v1 = volume_expanded(node1->bbox, bbox);
    if (v0 < v1) {
        *ibestnode = 0;
    } else {
        *ibestnode = 1;
    }
    return 0;
}

// volume_expanded returns the volume that bbox0 would be expanded by if it
// were merged with bbox1
float volume_expanded(Bbox *bbox0, Bbox *bbox1) {
    float v_original = (bbox0->xmax - bbox0->xmin) * (bbox0->ymax - bbox0->ymin);
    float v_new = (maxf(bbox0->xmax, bbox1->xmax) -
                    minf(bbox0->xmin, bbox1->xmin)) *
                  (maxf(bbox0->ymax, bbox1->ymax) -
                    minf(bbox0->ymin, bbox1->ymin));
    return v_new - v_original;
}

int is_within(Bbox *outerbb, Bbox *testbb) {
    if ((testbb->xmin > outerbb->xmin) && (testbb->xmax <= outerbb->xmax) &&
        (testbb->ymin > outerbb->ymin) && (testbb->ymax <= outerbb->ymax)) {
        return 1;
    } else {
        return 0;
    }
}

float bbox_intersection_area(Bbox *bb0, Bbox *bb1) {
    float dx = 0, dy = 0;
    dx = maxf(minf(bb0->xmax, bb1->xmax) - maxf(bb0->xmin, bb1->xmin), 0.0);
    dy = maxf(minf(bb0->ymax, bb1->ymax) - maxf(bb0->ymin, bb1->ymin), 0.0);
    return dx * dy;
}

int is_overlapping(Bbox *bb, Bbox *testbb) {
    if (bbox_intersection_area(bb, testbb) > 0.0) {
        return 1;
    } else {
        return 0;
    }
}
