""" Cython wrapper for rtree """

cdef extern from "rtree.h":

    cdef enum Strategy:
        LINEAR, QUADRATIC

    cdef enum NodeType:
        LEAF, NONLEAF

    cdef struct BoundingBox:
        float xmin
        float ymin
        float xmax
        float ymax
    ctypedef BoundingBox Bbox

    cdef struct RTreeNode:
        int count
        Bbox* bbox
    ctypedef RTreeNode Node

    cdef struct PointerPool:
        int size
        int count
    ctypedef PointerPool Pool

    Node* rt_new_node(NodeType, Strategy, int, Node*)
    Node* rt_insert(Node*, Bbox*, int)
    Bbox* rt_new_bbox()
    void rt_free(Node*)
    void print_bbox(Bbox*)
    Pool *rt_search_within(Node*, Bbox*, int)
    Pool *rt_search_overlapping(Node*, Bbox*, int)

    char *pool_pop(Pool*, int)
    void pool_destroy(Pool*)

cdef class RTree:
    cdef int count
    cdef Node* root

    def __init__(self, list geometries, maxchildren=50):
        self.root = NULL
        cdef int i = 0
        cdef object geom
        cdef Node* root

        root = rt_new_node(LEAF, LINEAR, maxchildren, NULL)
        for geom in geometries:
            if not hasattr(geom, "bbox"):
                raise AttributeError("cannot construct R-tree index from items "
                                     "missing a `bbox` attribute")
            _bb = geom.bbox
            bb = rt_new_bbox()
            bb.xmin = _bb[0]
            bb.ymin = _bb[1]
            bb.xmax = _bb[2]
            bb.ymax = _bb[3]
            root = rt_insert(root, bb, i)
            i += 1

        self.count = i
        self.root = root
        return

    def __dealloc__(self):
        if self.root != NULL:
            rt_free(self.root)

    @property
    def bbox(self):
        return (self.root.bbox.xmin, self.root.bbox.ymin,
                self.root.bbox.xmax, self.root.bbox.ymax)

    def search_within(self, bbox, int max_results=-1):
        """ Return a list of geometries that are within a bounding box. """
        cdef Pool* result
        cdef Bbox* bb = rt_new_bbox()
        cdef list out = []

        bb.xmin = bbox[0]
        bb.ymin = bbox[1]
        bb.xmax = bbox[2]
        bb.ymax = bbox[3]

        result = rt_search_within(self.root, bb, max_results)
        while result.count != 0:
            out.append((<int*> pool_pop(result, result.count-1))[0])
        pool_destroy(result)
        return out

    def search_overlapping(self, bbox, int max_results=-1):
        """ Return a list of geometries that are within a bounding box. """
        cdef Pool* result
        cdef Bbox* bb = rt_new_bbox()
        cdef list out = []

        bb.xmin = bbox[0]
        bb.ymin = bbox[1]
        bb.xmax = bbox[2]
        bb.ymax = bbox[3]

        result = rt_search_overlapping(self.root, bb, max_results)
        while result.count != 0:
            out.append((<int*> pool_pop(result, result.count-1))[0])
        pool_destroy(result)
        return out

