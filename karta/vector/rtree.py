"""
Implements an R-tree data structure
"""

import numpy as np

class Node(object):
    def __init__(self, bbox, parent, max_children=50):
        self.bbox = list(bbox)
        self.parent = parent
        self.max_children = max_children
        self.children = []

    def __len__(self):
        return len(self.children)

    def add(self, geom):
        raise NotImplementedError()

    def split(self):
        raise NotImplementedError()

class LeafNode(Node):
    def __init__(self, bbox, parent, **kw):
        super(LeafNode, self).__init__(bbox, parent, **kw)

    def add(self, geom):
        """ Add a geometry and return the possibly modified bbox """
        retval = 0
        newnodes = []
        bb = geom.bbox
        self.bbox = [min(self.bbox[0], bb[0]), min(self.bbox[1], bb[1]),
                     max(self.bbox[2], bb[2]), max(self.bbox[3], bb[3])]
        self.children.append(geom)

        if len(self.children) > self.max_children:
            retval = 1
            newnodes = self.split()
        return retval, newnodes

    def split(self, kind="linear"):
        """ Split in place and return a new node. """
        bbs = [c.bbox for c in self.children]
        seeds = linear_pick_seeds(self.children)
        node0 = LeafNode(bbs[seeds[0]], self.parent, max_children=self.max_children)
        node1 = LeafNode(bbs[seeds[1]], self.parent, max_children=self.max_children)

        # add each child to one of the two nodes
        children0 = []
        children1 = []

        for i, (bb, child) in enumerate(zip(bbs, self.children)):
            if volume_expanded(node0, child) < volume_expanded(node1, child):
                children0.append(child)
            else:
                children1.append(child)

        for child in children0:
            bb = child.bbox
            node0.bbox = [min(node0.bbox[0], bb[0]), min(node0.bbox[1], bb[1]), 
                          max(node0.bbox[2], bb[2]), max(node0.bbox[3], bb[3])]
            node0.children.append(child)

        for child in children1:
            bb = child.bbox
            node1.bbox = [min(node1.bbox[0], bb[0]), min(node1.bbox[1], bb[1]), 
                          max(node1.bbox[2], bb[2]), max(node1.bbox[3], bb[3])]
            node1.children.append(child)

        return [node0, node1]

class NonLeafNode(Node):
    def __init__(self, bbox, parent, **kw):
        super(NonLeafNode, self).__init__(bbox, parent, **kw)

    def add(self, geom):
        """ Add a geometry and return the possibly modified bbox """
        retval = 0
        newnodes = []

        # choose a child node to add to
        v = [volume_expanded(c, geom) for c in self.children]
        i = v.index(min(v))
        target = self.children[i]
        flag, newchildnodes = target.add(geom)
        bb = target.bbox
        self.bbox = [min(self.bbox[0], bb[0]), min(self.bbox[1], bb[1]),
                     max(self.bbox[2], bb[2]), max(self.bbox[3], bb[3])]

        if flag == 1:
            # replace
            del self.children[i]
            self.children.extend(newchildnodes)
            if len(self.children) > self.max_children:
                newnodes = self.split()
                retval = 1

        return retval, newnodes

    def split(self):
        bbs = [c.bbox for c in self.children]
        seeds = linear_pick_seeds(self.children)
        node0 = NonLeafNode(bbs[seeds[0]], self.parent, max_children=self.max_children)
        node1 = NonLeafNode(bbs[seeds[1]], self.parent, max_children=self.max_children)

        # add each child to one of the two nodes
        children0 = []
        children1 = []

        for i, (bb, child) in enumerate(zip(bbs, self.children)):
            if volume_expanded(node0, child) < volume_expanded(node1, child):
                children0.append(child)
            else:
                children1.append(child)

        for child in children0:
            bb = child.bbox
            node0.bbox = [min(node0.bbox[0], bb[0]), min(node0.bbox[1], bb[1]), 
                          max(node0.bbox[2], bb[2]), max(node0.bbox[3], bb[3])]
            node0.children.append(child)

        for child in children1:
            bb = child.bbox
            node1.bbox = [min(node1.bbox[0], bb[0]), min(node1.bbox[1], bb[1]), 
                          max(node1.bbox[2], bb[2]), max(node1.bbox[3], bb[3])]
            node1.children.append(child)

        return [node0, node1]

def build_tree(geoms, max_children=50):
    bboxes = [g.bbox for g in geoms]
    xmin = min(b[0] for b in bboxes)
    xmax = max(b[2] for b in bboxes)
    ymin = min(b[1] for b in bboxes)
    ymax = max(b[3] for b in bboxes)

    root = LeafNode([xmin, ymin, xmax, ymax], None, max_children=max_children)
    for geom in geoms:
        flag, newnodes = root.add(geom)
        if flag == 1:
            root = NonLeafNode(root.bbox, None, max_children=max_children)
            root.children = newnodes
    return root

def linear_pick_seeds(geoms):
    """ Choose the seeds for splitting a list of geometries using the linear
    approximation of Guttman (1984). """
    fullbb = list(geoms[0].bbox)
    extreme_x0 = 0  # highest xmin
    extreme_x1 = 0  # lowest xmax
    extreme_y0 = 0  # highest ymin
    extreme_y1 = 0  # lowest ymax
    i = 1
    for geom in geoms[1:]:
        gbb = geom.bbox
        fullbb[0] = min(fullbb[0], gbb[0])
        fullbb[1] = min(fullbb[1], gbb[1])
        fullbb[2] = max(fullbb[2], gbb[2])
        fullbb[3] = max(fullbb[3], gbb[3])
        if geoms[extreme_x0].bbox[0] < gbb[0]:
            extreme_x0 = i
        if geoms[extreme_x1].bbox[2] > gbb[2]:
            extreme_x1 = i
        if geoms[extreme_y0].bbox[1] < gbb[1]:
            extreme_y0 = i
        if geoms[extreme_y1].bbox[3] > gbb[3]:
            extreme_y1 = i
        i += 1
    # normalize
    dx = (geoms[extreme_x0].bbox[0] - geoms[extreme_x1].bbox[2]) / (fullbb[2]-fullbb[0])
    dy = (geoms[extreme_y0].bbox[1] - geoms[extreme_y1].bbox[3]) / (fullbb[3]-fullbb[1])
    if dx > dy:
        return extreme_x0, extreme_x1
    else:
        return extreme_y0, extreme_y1

def area(bbox):
    return (bbox[2]-bbox[0])*(bbox[3]-bbox[1])

def volume_expanded(geom0, geom1):
    """ Return the volume by which the union bounding box of geom0, geom1 is
    larger than the bounding box of geom0. """
    bb0 = geom0.bbox
    bb1 = geom1.bbox
    union_bbox = [min(bb0[0], bb1[0]), min(bb0[1], bb1[1]),
                  max(bb0[2], bb1[2]), max(bb0[3], bb1[3])]
    return area(union_bbox)-area(bb0)

def depth(tree, d=0):
    if isinstance(tree, LeafNode):
        return d+1
    else:
        if len(tree.children) != 0:
            return max(depth(child, d+1) for child in tree.children)
        else:
            return d

def itergeoms(tree):
    if isinstance(tree, LeafNode):
        for geom in tree.children:
            yield geom
    else:
        for child in tree.children:
            for geom in itergeoms(child):
                yield geom

def iterchildren(tree):
    for child in tree.children:
        if isinstance(child, Node):
            for _child in iterchildren(child):
                yield _child
        else:
            yield child

def iterleaves(tree):
    for child in tree.children:
        if isinstance(child, NonLeafNode):
            for _child in iterleaves(child):
                yield _child
    if isinstance(tree, LeafNode):
        yield tree
