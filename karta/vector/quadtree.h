#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pool.h"

typedef enum {LEAF, NONLEAF} NodeType;

typedef enum {UPPERLEFT, UPPERRIGHT, LOWERLEFT, LOWERRIGHT} Quadrant;

typedef struct BoundingBox {
    double xmin;
    double ymin;
    double xmax;
    double ymax;
} Bbox;

typedef struct PositionStruct {
    int id;
    double x;
    double y;
} Position;

typedef union NodePtrUnionTag {
    struct QuadTreeLeafNode *leafnode;
    struct QuadTreeNonleafNode *nonleafnode;
} NodePtrUnion;

typedef struct QuadTreeLeafNode {
    NodeType type;
    Bbox *bbox;
    int count;
    int max_positions;
    Position *positions;
    struct QuadTreeNonleafNode *parent;
} LeafNode;

typedef struct QuadTreeNonleafNode {
    NodeType type;
    Bbox *bbox;
    NodePtrUnion ulnode;
    NodePtrUnion urnode;
    NodePtrUnion llnode;
    NodePtrUnion lrnode;
    struct QuadTreeNonleafNode *parent;
} NonleafNode;

// forward declarations
void qt_free_node(NodePtrUnion);
NonleafNode *qt_split(LeafNode*);
Quadrant get_quadrant(Bbox*, Position*);
int isduplicate(LeafNode*, Position*);
int iswithin(Bbox*, Position*);
int overlaps(Bbox*, Bbox*);

Position* qt_new_position(int id, double x, double y) {
    Position *p = (Position*) malloc(sizeof(Position));
    p->id = id;
    p->x = x;
    p->y = y;
    return p;
}

LeafNode* qt_new_leaf(int max_positions, Bbox *bbox) {
    LeafNode *node = malloc(sizeof(LeafNode));
    node->type = LEAF;
    node->bbox = bbox;
    node->count = 0;
    node->max_positions = max_positions;
    node->positions = (Position*) malloc((max_positions+1) * sizeof(Position));
    node->parent = NULL;
    return node;
}

NonleafNode* qt_new_nonleaf(Bbox *bbox) {
    NonleafNode *node = malloc(sizeof(NonleafNode));
    node->type = NONLEAF;
    node->bbox = bbox;
    node->parent = NULL;
    return node;
}

Bbox* qt_new_bbox(double xmin, double ymin, double xmax, double ymax) {
    Bbox *bb = (Bbox*) malloc(sizeof(Bbox));
    bb->xmin = xmin;
    bb->ymin = ymin;
    bb->xmax = xmax;
    bb->ymax = ymax;
    return bb;
}

// Insert a position into a node. *dup* is a pointer to an integer that is set to
// a value in the range [0, max_positions) if a position is duplicated.
// Duplicate positions are not inserted, so the caller must keep track.
// returns NULL or a pointer to a new root node
NonleafNode *qt_insert(NodePtrUnion node_union, Position position, int *dup) {
    NonleafNode *retnode = NULL;
    NonleafNode *nonleaf = NULL;
    LeafNode *leaf = NULL;
    Quadrant quad;
    NodeType type = node_union.leafnode->type;
    int loops = 0;

    switch (type) {
        case LEAF: leaf = node_union.leafnode; break;
        case NONLEAF: nonleaf = node_union.nonleafnode; break;
    }

    while (1) {
        if (type == LEAF) {

            // check for duplicates - if so, break and return
            *dup = isduplicate(leaf, &position);
            if (*dup != -1) {
                break;
            }

            leaf->positions[leaf->count] = position;
            leaf->count++;
            if (leaf->count >= leaf->max_positions) {
                retnode = qt_split(leaf);
                if (loops != 0) {
                    retnode = NULL;
                }
            }
            break;
        }
        else if (type == NONLEAF) {

            quad = get_quadrant(nonleaf->bbox, &position);
            switch (quad) {
                case LOWERLEFT:
                    type = nonleaf->llnode.leafnode->type;
                    if (type == LEAF) {
                        leaf = nonleaf->llnode.leafnode;
                    }
                    else if (type == NONLEAF) {
                        nonleaf = nonleaf->llnode.nonleafnode;
                    }
                    break;

                case UPPERLEFT:
                    type = nonleaf->ulnode.leafnode->type;
                    if (type == LEAF) {
                        leaf = nonleaf->ulnode.leafnode;
                    }
                    else if (type == NONLEAF) {
                        nonleaf = nonleaf->ulnode.nonleafnode;
                    }
                    break;

                case LOWERRIGHT:
                    type = nonleaf->lrnode.leafnode->type;
                    if (type == LEAF) {
                        leaf = nonleaf->lrnode.leafnode;
                    }
                    else if (type == NONLEAF) {
                        nonleaf = nonleaf->lrnode.nonleafnode;
                    }
                    break;

                case UPPERRIGHT:
                    type = nonleaf->urnode.leafnode->type;
                    if (type == LEAF) {
                        leaf = nonleaf->urnode.leafnode;
                    }
                    else if (type == NONLEAF) {
                        nonleaf = nonleaf->urnode.nonleafnode;
                    }
                    break;
            }
        }
        else {
            printf("illegal node type\n");
            exit(1);
        }
        loops++;
    }
    return retnode;
}

// split a leaf node and return a reference to a new nonleaf node
NonleafNode *qt_split(LeafNode *node) {
    Bbox *bb = node->bbox;
    double xmid = 0.5*(bb->xmin + bb->xmax);
    double ymid = 0.5*(bb->ymin + bb->ymax);

    // create four new bboxes
    Bbox *bb_ul = qt_new_bbox(bb->xmin, ymid, xmid, bb->ymax);
    Bbox *bb_ur = qt_new_bbox(xmid, ymid, bb->xmax, bb->ymax);
    Bbox *bb_ll = qt_new_bbox(bb->xmin, bb->ymin, xmid, ymid);
    Bbox *bb_lr = qt_new_bbox(xmid, bb->ymin, bb->xmax, ymid);

    // create a new nonleaf parent
    NonleafNode *parent = qt_new_nonleaf(bb);

    // create new leaf nodes and modify input node
    LeafNode *node_ul = qt_new_leaf(node->max_positions, bb_ul);
    LeafNode *node_ur = qt_new_leaf(node->max_positions, bb_ur);
    LeafNode *node_ll = qt_new_leaf(node->max_positions, bb_ll);
    LeafNode *node_lr = qt_new_leaf(node->max_positions, bb_lr);

    // create relationships
    node_ul->parent = parent;
    node_ur->parent = parent;
    node_ll->parent = parent;
    node_lr->parent = parent;
    parent->ulnode.leafnode = node_ul;
    parent->urnode.leafnode = node_ur;
    parent->llnode.leafnode = node_ll;
    parent->lrnode.leafnode = node_lr;

    // connect parent to grandparent
    NonleafNode *grandparent = node->parent;
    Quadrant quad;
    Position p;
    if (grandparent != NULL) {
        p.x = xmid;
        p.y = ymid;
        quad = get_quadrant(grandparent->bbox, &p);
        switch (quad) {
            case LOWERLEFT: grandparent->llnode.nonleafnode = parent; break;
            case UPPERLEFT: grandparent->ulnode.nonleafnode = parent; break;
            case LOWERRIGHT: grandparent->lrnode.nonleafnode = parent; break;
            case UPPERRIGHT: grandparent->urnode.nonleafnode = parent; break;
        }
    }

    // distribute positions
    Position *pos;
    int i;
    for (i=0; i!=node->count; i++) {
        pos = &(node->positions[i]);
        quad = get_quadrant(bb, pos);

        switch (quad) {
            case LOWERLEFT:
                node_ll->positions[node_ll->count] = *pos;
                node_ll->count++;
                break;
            case LOWERRIGHT:
                node_lr->positions[node_lr->count] = *pos;
                node_lr->count++;
                break;
            case UPPERLEFT:
                node_ul->positions[node_ul->count] = *pos;
                node_ul->count++;
                break;
            case UPPERRIGHT:
                node_ur->positions[node_ur->count] = *pos;
                node_ur->count++;
                break;
        }
    }

    // check whether any child leaf is still over capacity
    if (node_ul->count >= node_ul->max_positions) {
        parent->ulnode.nonleafnode = qt_split(node_ul);
    }
    else if (node_ur->count >= node_ur->max_positions) {
        parent->urnode.nonleafnode = qt_split(node_ur);
    }
    else if (node_ll->count >= node_ll->max_positions) {
        parent->llnode.nonleafnode = qt_split(node_ll);
    }
    else if (node_lr->count >= node_lr->max_positions) {
        parent->lrnode.nonleafnode = qt_split(node_lr);
    }

    // destroy the old node, removing the bbox reference to preserve it
    node->bbox = NULL;
    NodePtrUnion node_union;
    node_union.leafnode = node;
    qt_free_node(node_union);
    return parent;
}

Pool *qt_search_within(NodePtrUnion node_union, Bbox *bbox) {
    Pool *results = pool_new(sizeof(int*), 64);
    Pool *quadrants = pool_new(sizeof(NodePtrUnion*), 16);
    NodePtrUnion active_node;
    NonleafNode *nonleaf;
    NodeType type;
    int i;

    pool_add(quadrants, (char*) &node_union);
    while (quadrants->count != 0) {

        active_node = *((NodePtrUnion*) pool_pop(quadrants, quadrants->count-1));
        type = active_node.leafnode->type;

        switch (type) {
            case LEAF:
                for (i=0; i!=active_node.leafnode->count; i++) {
                    if (iswithin(bbox, &(active_node.leafnode->positions[i])) == 1) {
                        pool_add(results, (char*) &(active_node.leafnode->positions[i]).id);
                    }
                }
                break;
            case NONLEAF:
                nonleaf = active_node.nonleafnode;
                if (overlaps(bbox, nonleaf->ulnode.leafnode->bbox)) {
                    pool_add(quadrants, (char*) &(nonleaf->ulnode));
                }
                if (overlaps(bbox, nonleaf->urnode.leafnode->bbox)) {
                    pool_add(quadrants, (char*) &(nonleaf->urnode));
                }
                if (overlaps(bbox, nonleaf->llnode.leafnode->bbox)) {
                    pool_add(quadrants, (char*) &(nonleaf->llnode));
                }
                if (overlaps(bbox, nonleaf->lrnode.leafnode->bbox)) {
                    pool_add(quadrants, (char*) &(nonleaf->lrnode));
                }
                break;
        }
    }

    pool_destroy(quadrants);
    return results;
}

void qt_free_position(Position *pos) {
    free(pos);
}

void qt_free_bbox(Bbox *bbox) {
    free(bbox);
}

void qt_free_node(NodePtrUnion node_union) {
    LeafNode *leafnode;
    NonleafNode *nonleafnode;
    NodeType type = node_union.leafnode->type;

    if (type == LEAF) {
        leafnode = node_union.leafnode;
        free(leafnode->bbox);
        free(leafnode->positions);
        free(leafnode);
    }
    else if (type == NONLEAF) {
        nonleafnode = node_union.nonleafnode;
        qt_free_node(nonleafnode->ulnode);
        qt_free_node(nonleafnode->urnode);
        qt_free_node(nonleafnode->llnode);
        qt_free_node(nonleafnode->lrnode);
        free(nonleafnode->bbox);
        free(nonleafnode);
    }
}

Quadrant get_quadrant(Bbox *bbox, Position *position) {
    double xmid = 0.5*(bbox->xmin + bbox->xmax);
    double ymid = 0.5*(bbox->ymin + bbox->ymax);
    if (position->x < xmid) {
        if (position->y < ymid) {
            return LOWERLEFT;
        } else {
            return UPPERLEFT;
        }
    } else {
        if (position->y < ymid) {
            return LOWERRIGHT;
        } else {
            return UPPERRIGHT;
        }
    }
}

int isduplicate(LeafNode *node, Position *pos) {
    for (int i=0; i!=node->count; i++) {
        if ((pos->x == node->positions[i].x) &&
            (pos->y == node->positions[i].y)) {
            return node->positions[i].id;
        }
    }
    return -1;
}

int iswithin(Bbox *bbox, Position *position) {
    if ((position->x > bbox->xmin) && (position->x < bbox->xmax) &&
        (position->y > bbox->ymin) && (position->y < bbox->ymax)) {
        return 1;
    }
    else {
        return 0;
    }
}

int overlaps(Bbox *bbox1, Bbox *bbox2) {
    if ((fmin(bbox1->xmax, bbox2->xmax) >= fmax(bbox1->xmin, bbox2->xmin)) &&
        (fmin(bbox1->ymax, bbox2->ymax) >= fmax(bbox1->ymin, bbox2->ymin))) {
        return 1;
    }
    else {
        return 0;
    }
}

