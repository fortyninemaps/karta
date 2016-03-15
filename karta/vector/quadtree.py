""" Implements a simple QuadTree datastructure. """

from ._cvectorgeo import iswithin, hashpt, bbox_intersection_area

class Node(object):
    def __init__(self, children, bbox, leaf):
        self.children = children
        self.bbox = bbox
        self.leaf = leaf

class QuadTree(object):
    """ Implements a convenience class that wraps a quadtree data structure and
    the methods needed to grow or query it.

    Parameters
    ----------
    bbox : 4-tuple
        (xmin, xmax, ymin, ymax).
    maxchildren : int, optional
        maximum number of points before a node is split
    maxdepth : int, optional
        maximum tree depth
    """

    def __init__(self, bbox, maxchildren=20, maxdepth=999):

        self.maxchildren = maxchildren
        self.maxdepth = maxdepth
        self.node = Node([], bbox, True)
        self.size = 0

    def addpt(self, pt):
        """ Add point to the QuadTree. Returns the depth at which the point was
        added. """
        (self.node, d) = addpt(self.node, pt, 1, self.maxchildren, self.maxdepth)
        self.size += 1
        return d

    def querypt(self, pt, method="hash"):
        """ Test whether QuadTree contains a point. *method* may be "recursion"
        [default] or "hash". """
        if method == "recursion":
            return querypt_recursion(self.node, pt)
        elif method == "hash":
            return querypt_hash(self.node, pt)
        else:
            raise ValueError("Unrecognized method '{0}'".format(method))

    def getfrombbox(self, bbox):
        """ Extract all points from a boundary box. """
        return getfrombbox(self.node, bbox)


def addpt(node, pt, depth, maxchildren, maxdepth):
    """ Add point to a node.

    Parameters
    ----------
    node : Node
    pt : tuple
    depth : int
        current quadtree level
    maxchildren : int
    maxdepth : int

    Returns
    -------
    Node
        updated node (which may or may not be the same as the starting node,
        depending on whether a split occurred)
    int
        depth at which point was inserted
    """
    if not iswithin(node.bbox, pt):
        raise BBoxError("({x},{y}) is not within bounds {bbox}".format(x=pt[0],
                                                                       y=pt[1],
                                                                       bbox=node.bbox))

    if node.leaf:

        if len(node.children) == maxchildren and depth != maxdepth:
            node = split(node)
            (node, depth) = addpt(node, pt, depth, maxchildren, maxdepth)

        else:
            node.children.append(pt)

    else:

        for (i,child) in enumerate(node.children):
            if iswithin(child.bbox, pt):
                (child, depth) = addpt(child, pt, depth+1, maxchildren, maxdepth)
                node.children[i] = child
                break

    return node, depth

def mean(x):
    n = len(x)
    if n != 0:
        return sum(x) / n
    else:
        None

def overlaps(bb0, bb1):
    return bbox_intersection_area(bb0, bb1) != 0.0

def split(node):
    """ Return a new node with the bbox of *node* and *node*'s children split
    between four child nodes. """
    (x0, y0, x1, y1) = node.bbox
    xm = 0.5 * (x0+x1)
    ym = 0.5 * (y0+y1)

    branch = Node([], node.bbox, False)
    bboxes = ((x0, y0, xm, ym), (xm, y0, x1, ym),
              (x0, ym, xm, y1), (xm, ym, x1, y1))

    for bbox in bboxes:
        pts = [pt for pt in node.children if iswithin(bbox, pt)]
        child = Node(pts, bbox, True)
        branch.children.append(child)
    return branch

def getfrombbox(parent, bbox):
    """ Return all points in *parent* node that are in *bbox*. """
    if parent.leaf:
        return [pt for pt in parent.children if iswithin(bbox, pt)]
    else:
        pts = []
        for child in parent.children:
            if overlaps(bbox, child.bbox):
                pts.extend(getfrombbox(child, bbox))
        return pts

def querypt_recursion(parent, pt):
    """ Test whether a point exists at *pt* using recursion. """
    if parent.leaf:
        for childpt in parent.children:
            if pt == childpt:
                return True
    else:
        for child in parent.children:
            if iswithin(child.bbox, pt):
                return querypt_recursion(child, pt)
    return False

def querypt_hash(parent, pt):
    (a, b, c, d) = parent.bbox
    h = hashpt(a, b, c, d, pt[0], pt[1])
    for quad in h:
        if parent.leaf:
            return pt in parent.children
        else:
            parent = parent.children[quad]

class BBoxError(Exception):
    def __init__(self, message):
        self.message = message
    def __repr__(self):
        return self.message
