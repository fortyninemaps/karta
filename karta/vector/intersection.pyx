from libc.math cimport NAN, M_PI, sin, cos, fmin, fmax, isnan
from cpython cimport bool
from coordstring cimport CoordString
from vectorgeo cimport (Vector2, Vector3,
                        cross2, cross3,
                        azimuth_sph, cart2sph, sph2cart,
                        eulerpole, eulerpole_cart)
import heapq

cdef bool isbetween_inc(double a, double b, double c):
    return fmin(a, c) <= b <= fmax(a, c)

cdef bool isbetween_incl(double a, double b, double c):
    return fmin(a, c) <= b < fmax(a, c)

cdef bool isbetween_incr(double a, double b, double c):
    return fmin(a, c) < b <= fmax(a, c)

ctypedef bool (*isbetween_t)(double, double, double)

def bboxes_overlap(tuple bb0, tuple bb1):
    """ Return whether planar bboxes overlap.
    """
    cdef float dx = 0.0
    cdef float dy = 0.0
    dx = fmin(bb0[2], bb1[2]) - fmax(bb0[0], bb1[0])
    dy = fmin(bb0[3], bb1[3]) - fmax(bb0[1], bb1[1])
    if dx == 0.0:
        dx = 1.0
    elif dx < 0.0:
        dx = 0.0
    if dy ==  0.0:
        dy = 1.0
    elif dy < 0.0:
        dy = 0.0
    if dx*dy == 0.0:
        return False
    else:
        return True

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
            if not isnan(xi):
                intersections.append((xi, yi))
    return intersections

cdef class Event(object):
    """ An Event represents an event in a sweep-line algorithm that signals a
    change in the sweep-line datastructure. This event signals the entrance and
    exits of line segments and segment intersections.
    """
    cdef public tuple seg
    cdef public int kind
    cdef public int color
    cdef public tuple left
    cdef public tuple right

    def __init__(self, tuple seg, int kind, int color):
        self.seg = seg
        self.kind = kind
        self.color = color
        self.left = None
        self.right = None

    def __richcmp__(self, Event other, int op):
        if op == 0:
            return self.seg < other.seg
        elif op == 1:
            return self.seg <= other.seg
        elif op == 2:
            return self.seg == other.seg
        elif op == 3:
            return self.seg != other.seg
        elif op == 4:
            return self.seg > other.seg
        elif op == 5:
            return self.seg >= other.seg
        else:
            raise ValueError("no __richcmp__ operation for op == %d" % op)

class EventQueuePlaceholder(object):
    """ This is a placeholder for an event-queue class. """

    def __init__(self):
        self.queue = []

    def __len__(self):
        return len(self.queue)

    def insert(self, double x, Event event):
        heapq.heappush(self.queue, (x, event))
        return

    def takefirst(self):
        """ return a flag and the first item from the queue. if flag == 1, then
        next item in the queue happens at the same time. """
        cdef int flag = 0
        cdef double x
        cdef Event event
        x, event = heapq.heappop(self.queue)
        if (len(self.queue) != 0) and (x == self.queue[0][0]):
            flag = 1
        return flag, x, event

class SweepLinePlaceholder(object):
    """ This is a placeholder for a sweep-line class. """
    def __init__(self, spherical=False):
        self.queue = []
        self.spherical = spherical

    def __iter__(self):
        cdef tuple item
        for item in self.queue:
            yield item

    def _get_slope_height(self, tuple seg, double x):
        cdef double m = 0.0, h = 0.0
        cdef tuple interx
        if self.spherical:
            m = -(azimuth_sph(Vector2(seg[0], seg[1]), Vector2(seg[2], seg[3]))*180/M_PI - 90)
            # compute intersection between seg and a meridional geodesic passing through x
            interx = intersection_sph(seg[0], seg[2], x, x, seg[1], seg[3], -90, 90)
            h = interx[1]
        else:
            if seg[3] == seg[1]:
                return m, seg[1]
            m = (seg[3]-seg[1]) / (seg[2]-seg[0])
            h = m*(x-seg[0]) + seg[1]
        return m, h

    def insert(self, tuple item, double x):
        if len(self.queue) == 0:
            self.queue.append(item)
            return

        cdef double m, y, m_, y_
        cdef int i = 0
        m, y = self._get_slope_height(item, x)
        m_, y_ = self._get_slope_height(self.queue[i], x)
        while (y_ < y) or ((y_ == y) and (m_ < m)):
            i += 1
            if i == len(self.queue):
                break
            m_, y_ = self._get_slope_height(self.queue[i], x)
        self.queue.insert(i, item)
        return

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

    def remove(self, item):
        cdef int idx
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

    cdef int na = len(a), nb = len(b)
    if a.ring:
        na = na + 1
    if b.ring:
        nb = nb + 1

    # initialize event queue
    #
    # event kind may be:
    # 0 : left point
    # 1 : right point
    # 2 : bottom point on vertical segment
    # 3 : top point on vertical segment
    # 4 : intersection (unused)
    cdef int i = 0
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
                event_queue.insert(x0, Event(seg, 3, 0))
            else:
                event_queue.insert(x0, Event(seg, 3, 0))
                event_queue.insert(x0, Event(seg, 2, 0))
        elif x0 < x1:
            event_queue.insert(x0, Event(seg, 0, 0))
            event_queue.insert(x1, Event(seg, 1, 0))
        else:
            event_queue.insert(x1, Event(seg, 0, 0))
            event_queue.insert(x0, Event(seg, 1, 0))

    for i in range(nb-1):
        x0 = b.getX(i)
        x1 = b.getX(i+1)
        y0 = b.getY(i)
        y1 = b.getY(i+1)
        seg = (x0, y0, x1, y1, 1)
        if x0 == x1:
            if y0 < y1:
                event_queue.insert(x0, Event(seg, 2, 1))
                event_queue.insert(x0, Event(seg, 3, 1))
            else:
                event_queue.insert(x0, Event(seg, 3, 1))
                event_queue.insert(x0, Event(seg, 2, 1))
        elif x0 < x1:
            event_queue.insert(x0, Event(seg, 0, 1))
            event_queue.insert(x1, Event(seg, 1, 1))
        else:
            event_queue.insert(x1, Event(seg, 0, 1))
            event_queue.insert(x0, Event(seg, 1, 1))

    # begin sweep search
    cdef double x
    cdef double yA
    cdef int event_type, segcolor
    cdef object sweepline = SweepLinePlaceholder()
    cdef Event event = None, ev = None
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

def intersects_sph(CoordString a, CoordString b):
    """ Shamos-Hoey intersection detection algorithm with adaptations for
    spherical coordinates. """

   # corner cases:
   # - simultaneous events          handled
   # - parallel overlapping lines   handled
   # - triple intersections         handled(?)
   # - vertical lines               handled
   # - self-intersections           handled
   #
   # all of the above need tests to verify

    cdef int na = len(a), nb = len(b)
    if a.ring:
        na = na + 1
    if b.ring:
        nb = nb + 1

    # initialize event queue
    #
    # event kind may be:
    # 0 : left point
    # 1 : right point
    # 2 : bottom point on vertical segment
    # 3 : top point on vertical segment
    # 4 : intersection (unused)
    cdef int i = 0
    cdef object event_queue = EventQueuePlaceholder()
    cdef double x0, x1, y0, y1
    cdef tuple seg

    for i in range(na-1):
        x0 = (a.getX(i) + 180.0) % 360.0 - 180.0
        x1 = (a.getX(i+1) + 180.0) % 360.0 - 180.0
        y0 = a.getY(i)
        y1 = a.getY(i+1)
        seg = (x0, y0, x1, y1, 0)
        if x0 == x1:
            if y0 < y1:
                event_queue.insert(x0, Event(seg, 2, 0))
                event_queue.insert(x0, Event(seg, 3, 0))
            else:
                event_queue.insert(x0, Event(seg, 3, 0))
                event_queue.insert(x0, Event(seg, 2, 0))
        elif x0 < x1:
            event_queue.insert(x0, Event(seg, 0, 0))
            event_queue.insert(x1, Event(seg, 1, 0))
        else:
            event_queue.insert(x1, Event(seg, 0, 0))
            event_queue.insert(x0, Event(seg, 1, 0))

    for i in range(nb-1):
        x0 = (b.getX(i) + 180.0) % 360.0 - 180.0
        x1 = (b.getX(i+1) + 180.0) % 360.0 - 180.0
        y0 = b.getY(i)
        y1 = b.getY(i+1)
        seg = (x0, y0, x1, y1, 1)
        if x0 == x1:
            if y0 < y1:
                event_queue.insert(x0, Event(seg, 2, 1))
                event_queue.insert(x0, Event(seg, 3, 1))
            else:
                event_queue.insert(x0, Event(seg, 3, 1))
                event_queue.insert(x0, Event(seg, 2, 1))
        elif x0 < x1:
            event_queue.insert(x0, Event(seg, 0, 1))
            event_queue.insert(x1, Event(seg, 1, 1))
        else:
            event_queue.insert(x1, Event(seg, 0, 1))
            event_queue.insert(x0, Event(seg, 1, 1))

    # begin sweep search
    cdef double x
    cdef double yA
    cdef int event_type, segcolor
    cdef object sweepline = SweepLinePlaceholder(spherical=True)
    cdef Event event = None, ev = None
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

                    yA = intersection_meridian(segA[0], segA[2], segA[1], segA[3], x)
                    if (seg[1] <= yA) and (seg[3] >= yA):
                        return True

                    if (seg[3] <= yA) and (seg[1] >= yA):
                        return True

                # can skip kind == 3 because it's covered by kind == 2

            elif event.kind == 0:
                # starting point
                segA = sweepline.leftof(seg)
                segB = sweepline.rightof(seg)
                if segA is not None and _intersects_sph(segA[0], segA[2],
                                                        seg[0], seg[2],
                                                        segA[1], segA[3],
                                                        seg[1], seg[3]):
                    return True
                if segB is not None and _intersects_sph(segB[0], segB[2],
                                                        seg[0], seg[2],
                                                        segB[1], segB[3],
                                                        seg[1], seg[3]):
                    return True

            elif event.kind == 1:
                # ending point
                segA = event.left
                segB = event.right
                if segA is not None and segB is not None \
                        and _intersects_sph(segB[0], segB[2], seg[0], seg[2],
                                            segB[1], segB[3], seg[1], seg[3]):
                    return True

    return False

cdef inline bool _intersects(tuple seg1, tuple seg2):
    """ Test whether two segments have an intersection. """
    cdef tuple interx
    if seg1[4] == seg2[4]:
        # same geometry
        return False
    if iscollinear(seg1[0], seg1[2], seg2[0], seg2[2],
            seg1[1], seg1[3], seg2[1], seg2[3]):
        return True
    interx = intersection(seg1[0], seg1[2], seg2[0], seg2[2],
                          seg1[1], seg1[3], seg2[1], seg2[3])
    if isnan(interx[0]):
        return False
    return True

#cdef inline Vector2 _project_gnomonic(Vector2 lonlat, lon0, R):
#    cdef double x = R*tan(lonlat.x - lon0)
#    cdef double y = R*tan(lonlat.y)/cos(lonlat.x-lon0)
#    return Vector2(x, y)
#
#cdef inline Vector2 _invproject_gnomonic(Vector2 xy, lon0, R):
#    return

cdef inline bool overlaps(double a0, double a1, double b0, double b1):
    if (a0 <= b0 <= a1) or (a0 <= b1 <= a1) or (b0 <= a0 <= b1) or (b0 <= a1 <= b1):
        return True
    else:
        return False

cdef bool iscollinear(double x0, double x1, double x2, double x3,
                      double y0, double y1, double y2, double y3):
    cdef double rxs = cross2(Vector2(x1-x0, y1-y0), Vector2(x3-x2, y3-y2))
    cdef double rxt = -1.0
    if rxs == 0:
        rxt = cross2(Vector2(x1-x0, y1-y0), Vector2(x3-x0, y3-y0))
        if rxt == 0:
            if (x1-x0) != 0 and overlaps(x0, x1, x2, x3):
                return True
            elif overlaps(y0, y1, y2, y3):
                return True
    return False

cdef bool iscollinear_sph(double x0, double x1, double x2, double x3,
                          double y0, double y1, double y2, double y3):
    cdef Vector3 ep1 = eulerpole(Vector2(x0, y0), Vector2(x1, y1))
    cdef Vector3 ep2 = eulerpole(Vector2(x2, y2), Vector2(x3, y3))
    # TODO: check overlap?
    if ((ep1.x == ep2.x) and (ep1.y == ep2.y) and (ep1.z == ep2.z)) \
            or ((ep1.x == -ep2.x) and (ep1.y == -ep2.y) and (ep1.z == -ep2.z)):
        return True
    else:
        return False

cpdef bool _intersects_sph(double x0, double x1, double x2, double x3,
                           double y0, double y1, double y2, double y3):
    cdef Vector3 ep1 = eulerpole(Vector2(x0, y0), Vector2(x1, y1))
    cdef Vector3 ep2 = eulerpole(Vector2(x2, y2), Vector2(x3, y3))
    cdef Vector2 normal = cart2sph(cross3(ep1, ep2))
    cdef double antipodal_lon = (normal.x + 360) % 360 - 180.0
    if isbetween_inc(x0, normal.x, x1) and isbetween_inc(x2, normal.x, x3):
        return True
    elif isbetween_inc(x0, antipodal_lon, x1) and isbetween_inc(x2, antipodal_lon, x3):
        return True
    else:
        return False

cdef double intersection_meridian(double x0, double x1, double y0, double y1,
                                  double xmeridian):
    """ Return the latitude of intersection of a spherical geodesic across a
    meridian. """
    cdef Vector3 ep1 = eulerpole(Vector2(x0, y0), Vector2(x1, y1))
    cdef Vector3 ep2 = Vector3(sin(M_PI*xmeridian/180.0),
                               cos(M_PI*xmeridian/180.0),
                               0.0)
    cdef Vector2 normal = cart2sph(cross3(ep1, ep2))
    normal.x = (normal.x + 180.0) % 360.0 - 180.0
    normal.y = (normal.y + 90.0) % 180.0 - 90.0
    cdef antipodal_lon = (normal.x + 360) % 360 - 180.0
    if isbetween_inc(x0, normal.x, x1):
        return normal.y
    elif isbetween_inc(x0, antipodal_lon, x1):
        return -normal.y
    else:
        return NAN

cpdef tuple intersection_sph(double x0, double x1, double x2, double x3,
                             double y0, double y1, double y2, double y3):
    """ Compute the spherical intersection between two segments,
    ((x0,y0), (x1,y1)) and ((x2,y2), (x3,y3)).

    Returns (NaN, NaN) if no intersection exists.
    """
    cdef Vector3 ep1 = eulerpole(Vector2(x0, y0), Vector2(x1, y1))
    cdef Vector3 ep2 = eulerpole(Vector2(x2, y2), Vector2(x3, y3))
    cdef Vector2 normal = cart2sph(cross3(ep1, ep2))
    normal.x = (normal.x + 180.0) % 360.0 - 180.0
    normal.y = (normal.y + 90.0) % 180.0 - 90.0
    cdef antipodal_lon = (normal.x + 360) % 360 - 180.0
    cdef Vector2 antipode = Vector2(antipodal_lon, -normal.y)
    if isbetween_inc(x0, normal.x, x1) and isbetween_inc(x2, normal.x, x3):
        return normal.x, normal.y
    elif isbetween_inc(x0, antipode.x, x1) and isbetween_inc(x2, antipode.x, x3):
        return antipode.x, antipode.y
    else:
        return NAN, NAN

cpdef tuple intersection(double x0, double x1, double x2, double x3,
                         double y0, double y1, double y2, double y3):
    """ Compute the planar intersection between two segments,
    ((x0,y0), (x1,y1)) and ((x2,y2), (x3,y3)).

    Returns (NaN, NaN) if no intersection exists.
    """
    cdef double rxs = cross2(Vector2(x1-x0, y1-y0), Vector2(x3-x2, y3-y2))
    if rxs == 0:
        return NAN, NAN

    cdef double t = cross2(Vector2(x2-x0, y2-y0), Vector2(x3-x2, y3-y2)) / rxs
    cdef double u = cross2(Vector2(x2-x0, y2-y0), Vector2(x1-x0, y1-y0)) / rxs
    if (0 < t <= 1) and (0 < u <= 1):
        return x0 + t*(x1-x0), y0 + t*(y1-y0)
    else:
        return NAN, NAN

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

