import numpy as np
cimport numpy as np
from coordstring cimport CoordString
from cpython cimport bool

cdef inline double dbl_max(double a, double b): return a if a >= b else b
cdef inline double dbl_min(double a, double b): return a if a <= b else b

cdef bool isbetween_inc(double a, double b, double c):
    return dbl_min(a, c) <= b <= dbl_max(a, c)

cdef bool isbetween_incl(double a, double b, double c):
    return dbl_min(a, c) <= b < dbl_max(a, c)

cdef bool isbetween_incr(double a, double b, double c):
    return dbl_min(a, c) < b <= dbl_max(a, c)

ctypedef bool (*isbetween_t)(double, double, double)

cdef inline double cross2(double u0, double u1, double v0, double v1):
    return u0*v1 - u1*v0

def all_intersections(CoordString a, CoordString b):
    """ Brute-force intersection search with len(a)*len(b) complexity """
    cdef int na = len(a)
    cdef int nb = len(b)
    cdef int i = 0, j = 0
    cdef double x0, x1, x2, x3, y0, y1, y2, y3
    cdef xi, yi
    cdef list intersections = []

    if not a.ring:
        na -= 1
    if not b.ring:
        nb -= 1

    for i in range(na):
        x0 = a.getX(i)
        x1 = a.getX(i+1)
        y0 = a.getY(i)
        y1 = a.getY(i+1)
        for j in range(nb):
            x2 = b.getX(j)
            x3 = b.getX(j+1)
            y2 = b.getY(j)
            y3 = b.getY(j+1)
            xi, yi = intersection(x0, x1, x2, x3, y0, y1, y2, y3)
            if not np.isnan(xi):
                intersections.append((xi, yi))
    return intersections

class Event(object):
    def __init__(self, seg, kind, color):
        self.seg = seg
        self.kind = kind
        self.color = color
        self.left = None
        self.right = None
    def __lt__(self, other):
        return self.seg < other.seg

class EventQueuePlaceholder(object):
    """ This is a placeholder for an event-queue class. It's not efficient, but
    it works. """
    def __init__(self):
        self.queue = []

    def __len__(self):
        return len(self.queue)

    def insert(self, double x, object event):
        self.queue.append((x, event))
        self.queue.sort()
        return

    def takefirst(self):
        """ return a flag and the first item from the queue. if flag == 1, then
        next item in the queue happens at the same time. """
        cdef int flag = 0
        cdef double x
        cdef object event
        x, event = self.queue[0]
        self.queue = self.queue[1:]
        if (len(self.queue) != 0) and (x == self.queue[0][0]):
            flag = 1
        return flag, x, event

class SweepLinePlaceholder(object):
    """ This is a placeholder for a sweep-line class. """
    def __init__(self):
        self.queue = []

    def __iter__(self):
        cdef tuple item
        for item in self.queue:
            yield item

    def _get_height(self, tuple seg, double x):
        if seg[2] == seg[1]:
            return seg[1]
        cdef double m
        m = (seg[3]-seg[1]) / (seg[2]-seg[0])
        return m*(x-seg[0]) + seg[1]

    def leftof(self, tuple item):
        cdef int idx
        idx = self.queue.index(item)
        if idx == 0:
            return None
        else:
            return self.queue[idx-1]

    def rightof(self, tuple item):
        cdef int idx
        idx = self.queue.index(item)
        if idx == len(self.queue) - 1:
            return None
        else:
            return self.queue[idx+1]

    def insert(self, tuple item, double x):
        if len(self.queue) == 0:
            self.queue.append(item)
            return

        y = self._get_height(item, x)
        i = 0
        while self._get_height(self.queue[i], x) < y:
            i += 1
            if i == len(self.queue):
                break
        self.queue.insert(i, item)
        return

    def remove(self, item):
        idx = self.queue.index(item)
        self.queue.pop(idx)
        return

def intersects(CoordString a, CoordString b):
    """ Shamos-Hoey intersection detection algorithm """

   # corner cases:
   # - simultaneous events          handled
   # - parallel overlapping lines   handled
   # - triple intersections         handled(?)
   # - vertical lines               handled
   # - self-intersections           handled
   #
   # all of the above need tests to verify

    cdef int na, nb, i = 0
    na = len(a)
    nb = len(b)
    if a.ring:
        na += 1
    if b.ring:
        nb += 1

    # initialize event queue
    #
    # event kind may be:
    # 0 : left point
    # 1 : right point
    # 2 : bottom point on vertical segment
    # 3 : top point on vertical segment
    # 4 : intersection (unused)
    cdef object event_queue = EventQueuePlaceholder()
    cdef double x0, x1, y0, y1
    cdef tuple seg

    for i in range(na-1):
        x0 = a.getX(i)
        x1 = a.getX(i+1)
        y0 = a.getY(i)
        y1 = a.getY(i+1)
        seg = (x0, y0, x1, y1, 0)
        if x0 == x1:
            if y0 < y1:
                event_queue.insert(x0, Event(seg, 2, 0))
                event_queue.insert(x1, Event(seg, 3, 0))
            else:
                event_queue.insert(x0, Event(seg, 3, 0))
                event_queue.insert(x1, Event(seg, 2, 0))
        elif x0 < x1:
            event_queue.insert(x0, Event(seg, 0, 0))
            event_queue.insert(x1, Event(seg, 1, 0))
        else:
            event_queue.insert(x0, Event(seg, 1, 0))
            event_queue.insert(x1, Event(seg, 0, 0))

    for i in range(nb-1):
        x0 = b.getX(i)
        x1 = b.getX(i+1)
        y0 = b.getY(i)
        y1 = b.getY(i+1)
        seg = (x0, y0, x1, y1, 1)
        if x0 == x1:
            if y0 < y1:
                event_queue.insert(x0, Event(seg, 2, 1))
                event_queue.insert(x1, Event(seg, 3, 1))
            else:
                event_queue.insert(x0, Event(seg, 3, 1))
                event_queue.insert(x1, Event(seg, 2, 1))
        elif x0 < x1:
            event_queue.insert(x0, Event(seg, 0, 1))
            event_queue.insert(x1, Event(seg, 1, 1))
        else:
            event_queue.insert(x0, Event(seg, 1, 1))
            event_queue.insert(x1, Event(seg, 0, 1))

    # begin sweep search
    cdef double x
    cdef double yA
    cdef int event_type, segcolor
    cdef object sweepline = SweepLinePlaceholder()
    cdef object event = None, ev = None
    cdef tuple segA, segB
    cdef list current_events = [], neighbours = []

    while len(event_queue) != 0:

        # get all simultaneous events
        current_events = []
        flag = 1
        while flag == 1:
            flag, x, event = event_queue.takefirst()
            current_events.append((x, event))

        # for all segments that begin, insert into sweepline
        for x, event in current_events:
            if event.kind == 0:
                sweepline.insert(event.seg, x)

        # for all segments that end, get the neighbours that are not also
        # ending and remove from sweepline
        for _, event in current_events:
            if event.kind == 1:
                segA = event.seg
                segB = event.seg
                while segA is not None and \
                        segA in (ev.seg for x, ev in current_events if ev.kind==1):
                    segA = sweepline.leftof(segA)
                while segB is not None and \
                        segB in (ev.seg for x, ev in current_events if ev.kind==1):
                    segB = sweepline.rightof(segB)
                event.left = segA
                event.right = segB
                sweepline.remove(event.seg)

        # check for intersections due to last set of events
        for x, event in current_events:

            seg = event.seg

            if (event.kind == 2):
                # vertical segment - check whether span crosses any segments
                for segA in sweepline:

                    if seg[4] == segA[4]:
                        continue

                    yA = (segA[3]-segA[1]) / (segA[2]-segA[0])*(seg[0]-segA[0]) + segA[1]
                    if (seg[1] <= yA) and (seg[3] >= yA):
                        return True

                    if (seg[3] <= yA) and (seg[1] >= yA):
                        return True

                # can skip kind == 3 because it's covered by kind == 2

            elif event.kind == 0:
                # starting point
                segA = sweepline.leftof(seg)
                segB = sweepline.rightof(seg)
                if segA is not None and _intersects(segA, seg):
                    return True
                if segB is not None and _intersects(seg, segB):
                    return True

            elif event.kind == 1:
                # ending point
                segA = event.left
                segB = event.right
                if segA is not None and segB is not None and _intersects(segA, segB):
                    return True

    return False

cdef inline bool _intersects(tuple seg1, tuple seg2):
    cdef tuple interx
    if seg1[4] == seg2[4]:
        # same geometry
        return False
    if iscollinear(seg1[0], seg1[2], seg2[0], seg2[2],
            seg1[1], seg1[3], seg2[1], seg2[3]):
        return True
    interx = intersection(seg1[0], seg1[2], seg2[0], seg2[2],
                          seg1[1], seg1[3], seg2[1], seg2[3])
    if np.isnan(interx[0]):
        return False
    return True

# def all_intersections_(CoordString a, CoordString b):
#     """ Bentley-Ottman intersection sweep with (n+k) log n complexity, where n
#     = len(a)+len(b) and k is the number of intersections. """
# 
#     # There are various corner cases to handle:
#     #
#     # - multiple events at the same place in the event queue
#     # - parallel overlapping lines
#     # - triple intersections
#     # - vertical lines
#     # - segment vertices
# 
#     cdef list event_queue, output
#     cdef int na, nb, segnum
#     cdef int i
# 
#     na = len(a)
#     nb = len(b)
#     if a.ring:
#         a += 1
#     if b.ring:
#         b += 1
# 
#     # initialize event queue
#     # values are tuple of (x, y, LINE, SEGA, SEGB)
#     # if LINE == 0, event is on line A
#     # if LINE == 1, event is on line B
#     # if LINE == -1, event is an intersection between SEGA, SEGB
#     event_queue = EventQueue()
#     segnum = 0
#     for i in range(na-1):
#         event_queue.insert((a.getX(i), a.getY(i), 0, segnum, 0))
#         event_queue.insert((a.getX(i+1), a.getY(i+1), 0, segnum, 0))
#         segnum += 1
# 
#     for i in range(nb-1):
#         event_queue.insert((b.getX(i), b.getY(i), 1, segnum, 0))
#         event_queue.insert((b.getX(i+1), b.getY(i+1), 1, segnum, 0))
#         segnum += 1
# 
#     sweepline = SweepTree()
#     output = []
# 
#     cdef tuple seg1, seg2, event, interx
#     cdef int segAbove, segBelow
#     cdef Node N1, N2, NAbove, NBelow
# 
#     while len(event_queue) != 0:
# 
#         event = event.pop()
#         if event[2] != -1:
#             seg1 = segments[event[3]]
#             if (event[0] == seg1[0]) and (event[1] == seg1[1]):
#                 # entering point
#                 N1 = sweepline.insert(seg1, event[0])
#                 NAbove = N1.one_right
#                 NBelow = N1.one_left
#                 _intersection_to_queue(seg1, NAbove.seg, event_queue)
#                 _intersection_to_queue(seg1, NBelow.seg, event_queue)
# 
#             else:
#                 # leaving point
#                 N1 = sweepline.find(seg1, event[0])
#                 NAbove = N1.one_right
#                 NBelow = N1.one_left
#                 sweepline.delete(seg1, event[0])
#                 _intersection_to_queue(NAbove.seg, NBelow.seg, event_queue)
# 
#         else:
#             # intersection point
#             seg1 = segments[event[3]]
#             seg2 = segments[event[4]]
#             N1 = sweepline.find(seg1, event[0])
#             N2 = sweepline.find(seg2, event[0])
# 
#             # swap segment ordering
#             N1.segment, N2.segment = N2.segment, N1.segment
# 
#             # establish which segment is above the other
#             if slope(N1.segment) > slope(N2.segment):
#                 # N1 is above (right)
#                 NAbove = N1.one_right
#                 NBelow = N2.one_left
# 
#                 _intersection_to_queue(N1.seg, NAbove.seg, event_queue)
#                 _intersection_to_queue(N2.seg, NBelow.seg, event_queue)
# 
#             else:
#                 # N1 is below (left)
#                 NBelow = N1.one_left
#                 NAbove = N2.one_right
# 
#                 _intersection_to_queue(N1.seg, NBelow.seg, event_queue)
#                 _intersection_to_queue(N2.seg, NAbove.seg, event_queue)

# cdef _intersection_to_queue(tuple seg1, tuple seg2, RedBlackTree queue):
#     """ Check whether seg1 and seg2 intersect. If so, attempt to add them to a
#     queue, returning the result. If not, return -1. """
#     cdef tuple interx, event
#     interx = intersection(seg1[0], seg1[2], seg2[0], seg2[2],
#                           seg1[1], seg1[3], seg2[1], seg2[3])
#     if not np.isnan(interx[0]):
#         event = (interx[0], interx[1], -1, N1.index, NBelow.index)
#         return queue.insert(event)
#     return -1

cdef inline bool overlaps(double a0, double a1, double b0, double b1):
    if (a0 <= b0 <= a1) or (a0 <= b1 <= a1) or (b0 <= a0 <= b1) or (b0 <= a1 <= b1):
        return True
    else:
        return False

cdef bool iscollinear(double x0, double x1, double x2, double x3,
                      double y0, double y1, double y2, double y3):
    cdef double rxs = cross2(x1-x0, y1-y0, x3-x2, y3-y2)
    cdef double rxt = -1.0
    if rxs == 0:
        rxt = cross2(x1-x0, y1-y0, x3-x0, y3-y0)
        if rxt == 0:
            if (x1-x0) != 0 and overlaps(x0, x1, x2, x3):
                return True
            elif overlaps(y0, y1, y2, y3):
                return True
    return False

cpdef intersection(double x0, double x1, double x2, double x3,
                   double y0, double y1, double y2, double y3):
    """ Returns coordinates of intersection point, or (NaN, NaN) if lines don't
    intersect.

    Line1 consists of point pairs 0, 1
    Line2 consists of point pairs 2, 3
    """
    cdef double rxs = cross2(x1-x0, y1-y0, x3-x2, y3-y2)
    if rxs == 0:
        # parallel or collinear
        return np.nan, np.nan

    cdef double t = cross2(x2-x0, y2-y0, x3-x2, y3-y2) / rxs
    cdef double u = cross2(x2-x0, y2-y0, x1-x0, y1-y0) / rxs
    if (0 < t <= 1) and (0 < u <= 1):
        return x0 + t*(x1-x0), y0 + t*(y1-y0)
    else:
        # non-intersecting
        return np.nan, np.nan

def count_crossings(double xp, double yp, CoordString coords):
    """ Count the number of times a vertical ray from (xp, yp) crosses a
    line/polygon defined by *coords*
    """
    cdef int n = len(coords)
    cdef int i, cnt = 0
    if not coords.ring:
        n -= 1
    cdef double x0 = coords.getX(0)
    cdef double y0 = coords.getY(0)
    cdef double x1, y1
    for i in range(1, n):
        x1 = coords.getX(i)
        y1 = coords.getY(i)
        if intersects_cn(xp, yp, x0, x1, y0, y1) == 1:
            cnt += 1
        x0 = x1
        y0 = y1
    return cnt

cdef int intersects_cn(double xp, double yp, double x0, double x1, double y0, double y1):
    """ Test whether a vertical ray emanating up from a point (xp, yp) crosses
    a line segment [(x0, y0), (x1, y1)]. Used to implement a crossing number
    membership test.
    """
    cdef double m, y

    if x0 != x1:
        m = (y1-y0) / (x1-x0)
    else:
        return 0

    y = y0 + m * (xp - x0)

    if y < yp:
        return 0

    cdef int iswithinx = 0
    cdef int iswithiny = 0

    if m > 0.0 and isbetween_incr(y0, y, y1):
        iswithiny = 1
    elif isbetween_incl(y0, y, y1):
        iswithiny = 1
    elif abs(y0-y1) < 1e-15 and abs(y-y0) < 1e-15:
        iswithiny = 1

    if isbetween_incr(x0, xp, x1):
        iswithinx = 1

    return iswithinx*iswithiny

