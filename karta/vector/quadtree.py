""" Implements a simple quadtree datastructure, with emphasis on performance. """

from collections import namedtuple

Node = namedtuple("Node", ["children", "bbox", "leaf"])

class QuadTree(object):
    """ Implements a convenience class that wraps a quadtree data structure and
    the methods needs to grow or query it. """

    def __init__(self, bbox, maxchildren=20, maxdepth=999):

        self.maxchildren = maxchildren
        self.maxdepth = maxdepth
        self.node = Node([], bbox, True)
        self.size = 0

    def addpt(self, pt):
        self.node = addpt(self.node, pt, 1, self.maxchildren, self.maxdepth)
        self.size += 1

    def querypt(self, pt):
        return querypt(self.node, pt)

    def getfrombbox(self, bbox):
        return getfrombbox(self.node, bbox)


def addpt(node, pt, depth, maxchildren, maxdepth):
    """ Add *pt* to *node*. *depth* is the current level in the quadtree.
    Returns the updated node (which may or may not be the same as the starting
    node, depending on whether a split occurred, as well as the depth at which
    the point was inserted. """

    if node.leaf:

        if len(node.children) == maxchildren and depth != maxdepth:
            node = split(node)
            addpt(node, pt, depth, maxchildren, maxdepth)

        else:
            node.children.append(pt)

    else:

        for child in node.children:
            if iswithin(child.bbox, pt):
                (node, depth) = addpt(child, pt, depth+1, maxchildren, maxdepth)
                break

    return node, depth

def mean(x):
    n = len(x)
    if n != 0:
        return sum(x) / n
    else:
        None

def iswithin(bbox, pt):
    if (bbox[0] <= pt[0] < bbox[1] and bbox[2] <= pt[1] < bbox[3]):
        return True
    else:
        return False

def overlaps(bbox0, bbox1):
    pts0 = ((x, y) for x in bbox0[:2] for y in bbox[2:])
    if True in (iswithin(bbox1, pt) for pt in pts0):
        return True
    pts1 = ((x, y) for x in bbox1[:2] for y in bbox[2:])
    if True in (iswithin(bbox0, pt) for pt in pts1):
        return True
    return False

def split(node):
    """ Return a new node with the bbox of *node* and *node*'s children split
    between four child nodes. """
    (x0, x1, y0, y1) = node.bbox
    xm = mean([pt[0] for pt in node.children])
    ym = mean([pt[1] for pt in node.children])

    branch = Node([], node.bbox, False)
    eps = 3e-15
    bboxes = ((x0, xm, y0, ym), (xm, x1+eps, y0, ym),
              (x0, xm, ym, y1+eps), (xm, x1+eps, ym, y1+eps))

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

def querypt(parent, pt):
    """ Test whether a point exists at *pt*. """
    if not parent.leaf:
        for child in parent.children:
            if iswithin(child.bbox, pt):
                return querypt(child, pt)
    else:
        for childpt in parent.children:
            if pt == childpt:
                return True
    return False

