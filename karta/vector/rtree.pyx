""" Cython wrapper for rtree """

cdef enum Strategy:
    LINEAR, QUADRATIC

cdef enum NodeType:
    LEAF, NONLEAAF

cdef struct BoundingBox:
    float xmin
    float ymin
    float xmax
    float ymax

ctypedef BoundingBox Bbox

cdef struct RTreeNode:
    NodeType type
    Strategy strategy
    int count
    int maxchildren
    RTreeNode* parent
    char** children
    int* indices
    Bbox* bbox

ctypedef RTreeNode Node

cdef extern from "rtree.h":
    Node* rt_new_node(NodeType, Strategy, int, Node*)
    Node* rt_insert(Node*, Bbox*, int)
    Bbox* rt_new_bbox()
    void rt_free(Node*)
    void print_bbox(Bbox*)

cdef class RTree:
    cdef list geometries
    cdef Node* root

    def __init__(self, list geometries, maxchildren=50):
        cdef int i
        cdef object geom
        cdef Node* root

        root = rt_new_node(LEAF, LINEAR, maxchildren, NULL)
        for i, geom in enumerate(geometries):
            _bb = geom.bbox
            bb = rt_new_bbox()
            bb.xmin = _bb[0]
            bb.ymin = _bb[1]
            bb.xmax = _bb[2]
            bb.ymax = _bb[3]
            root = rt_insert(root, bb, i)

        self.geometries = geometries
        self.root = root
        return

    def __dealloc__(self):
        rt_free(self.root)

    @property
    def bbox(self):
        return (self.root.bbox.xmin, self.root.bbox.ymin,
                self.root.bbox.xmax, self.root.bbox.ymax)

    def search(self):
        raise NotImplementedError()

