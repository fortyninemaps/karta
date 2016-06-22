""" Cython wrapper quadtree """

# from .coordstring import Coordstring
# from .coordstring cimport Coordstring

cdef extern from "quadtree.h":

    ctypedef enum NodeType:
        LEAF, NONLEAF

    ctypedef enum Quadrant:
        UPPERLEFT, UPPERRIGHT, LWERLEFT, LOWERRIGHT

    struct BoundingBox:
        double xmin
        double ymin
        double xmax
        double ymax
    ctypedef BoundingBox Bbox

    struct PositionStruct:
        int id
        double x
        double y
    ctypedef PositionStruct Position

    struct QuadTreeLeafNode:
        Bbox *bbox
    ctypedef QuadTreeLeafNode LeafNode

    struct QuadTreeNonleafNode:
        Bbox *bbox
    ctypedef QuadTreeNonleafNode NonleafNode

    union NodePtrUnionTag:
        QuadTreeLeafNode *leafnode
        QuadTreeNonleafNode *nonleafnode
    ctypedef NodePtrUnionTag NodePtrUnion

    cdef struct PointerPool:
        int size
        int count
    ctypedef PointerPool Pool

    Position *qt_new_position(int, double, double)
    void qt_free_position(Position*)

    Bbox *qt_new_bbox(double, double, double, double)
    void qt_free_bbox(Bbox*)

    LeafNode *qt_new_leaf(int, Bbox*)
    NonleafNode *qt_new_nonleaf(Bbox*)
    NonleafNode *qt_insert(NodePtrUnion, Position, int*)
    void qt_free_node(NodePtrUnion)
    Pool *qt_search_within(NodePtrUnion, Bbox*)
    char *pool_pop(Pool*, int)

cdef class QuadTree:
    cdef int count
    cdef NodePtrUnion root
    cdef readonly dict duplicates

    def __cinit__(self, points, int leaf_capacity=50):
        cdef NonleafNode *retnode
        cdef Position *pos
        cdef Bbox *bbox
        cdef double x, y
        cdef double xmin, ymin, xmax, ymax
        cdef int i = 0
        cdef int idup = 0
        cdef dict duplicates = {}

        xmin, ymin, xmax, ymax = points.bbox
        bbox = qt_new_bbox(xmin, ymin, xmax, ymax)
        self.root.leafnode = qt_new_leaf(leaf_capacity, bbox)

        while i != len(points):
            x, y = points[i][:2]
            pos = qt_new_position(i, x, y)
            retnode = qt_insert(self.root, pos[0], &idup)
            if idup != -1:
                if idup in duplicates:
                    duplicates[idup].append(i)
                else:
                    duplicates[idup] = [i]
            elif (retnode != NULL):
                self.root.nonleafnode = retnode
            qt_free_position(pos)
            i += 1

        self.count = i
        self.duplicates = duplicates
        return

    def __dealloc__(self):
        if not (self.root.leafnode == NULL):
            qt_free_node(self.root)

    def __len__(self):
        return self.count

    @property
    def bbox(self):
        cdef NonleafNode *node = self.root.nonleafnode
        cdef Bbox *bbox = node.bbox
        cdef double xmin = bbox.xmin
        cdef double ymin = bbox.ymin
        cdef double xmax = bbox.xmax
        cdef double ymax = bbox.ymax
        return (xmin, ymin, xmax, ymax)

    def search_within(self, double xmin, double ymin, double xmax, double ymax):
        cdef Bbox *bbox = qt_new_bbox(xmin, ymin, xmax, ymax)
        cdef Pool *results = qt_search_within(self.root, bbox)
        cdef int *idx
        cdef list out = []

        while results.count != 0:
            idx = <int*> pool_pop(results, results.count-1)
            if idx[0] in self.duplicates:
                out.extend(self.duplicates[idx[0]])
            out.append(idx[0])
        qt_free_bbox(bbox)
        return out

