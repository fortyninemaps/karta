def interpolate1(float x, float y, float a, float b, float c, float d):
    """ Return a value *v(x,y)* in the regular structured stencil
    
            a --- b
            |  v  |
            c --- d
    
    using linear interpolation. The coordinates (x, y) must be normalized by
    the horizontal and vertical grid spacings, respectively.
    """
    cdef float left, right, res
    
    left = (c-a)*y + a
    right = (d-b)*y + b
    res = (right - left) * x + left
    return res
