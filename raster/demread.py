""" Functions for reading the USGS DEM format. """

import itertools

def coerce_float(s):
    return float(s.replace("D","e",1)
                  .replace("E","e",1)
                  .replace("F","e",1)
                  .replace("G","e",1))

FMT_KEY = {"A":str,
           "I":int,
           "D":coerce_float,
           "E":coerce_float,
           "F":coerce_float,
           "G":coerce_float}


BLOCKA = [("fnm",0,40,"A40"),
          ("free",40,109,"A69"),
          ("se_corner",109,135,"2(I4,I2,F7.4)"),
          ("proc_code",135,136,"A1",),
          ("sectional_indicator",137,140,"A3"),
          ("origin_code",140,144,"A4"),
          ("dem_level_code",144,150,"I6"),
          ("elev_pattern",150,156,"I6"),
          ("crs",156,162,"I6"),
          ("crs_zone",162,168,"I6"),
          ("proj",168,528,"15D24.15"),
          ("crs_unit",528,534,"I6"),
          ("z_unit",534,540,"I6"),
          ("nsides",540,546,"I6"),
          ("quad_coords",546,738,"4(2D24.15)"),
          ("minmax",738,786,"2D24.15"),
          ("angle_ccw",786,810,"D24.15"),
          ("z_acc",810,816,"I6"),
          ("resolution",816,852,"3E12.6"),
          ("size",852,864,"2I6"),       # Old format end
          ("max_contour_int",864,869,"I5"),
          ("max_contour_units",869,870,"I1"),
          ("min_contour_int",870,875,"I5"),
          ("min_contour_units",875,876,"I1"),
          ("src_year",876,880,"I4"),
          ("rev_year",880,884,"I4"),
          ("inspection_flag",884,885,"A1"),
          ("validation_flag",885,886,"A1"),
          ("void_flag",886,888,"I1"),
          ("z_datum",888,890,"I2"),
          ("xy_datum",890,892,"I2"),
          ("edition",892,896,"I4"),
          ("void_perc",896,900,"I4"),
          ("edge_flag",900,908,"4I2"),
          ("z_shift",908,915,"F7.2")]

def nreps(fmt):
    """ Return the number of characters in a single record of
    a potentially multiple-record string. """
    try:
        reps = int(fmt[0])
        return reps * nreps(fmt[1:].lstrip("(").rstrip(")"))

    except ValueError:
        return 1
        #return fmt.count(",") + 1

def reclen(fmt):
    """ Return the number of characters in a single record of
    a potentially multiple-record string. """
    try:
        reps = int(fmt[0])
        return reclen(fmt[1:].lstrip("(").rstrip(")"))

    except ValueError:
        if "," in fmt:
            return map(lambda a: reclen(a)[0], fmt.split(","))

        if fmt[0] in ("G", "F", "E", "D"):
            nch = int(fmt[1:].split(".")[0])
        elif fmt[0] in ("A", "I"):
            nch = int(fmt[1:])
        return [nch]

def dtype(fmt):
    """ Return a constructor for the datatype encoded in *fmt*. """
    try:
        reps = int(fmt[0])
        return dtype(fmt[1:].lstrip("(").rstrip(")"))

    except ValueError:
        if "," in fmt:
            return map(lambda a: dtype(a)[0], fmt.split(","))
        else:
            return [FMT_KEY[fmt[0]]]

def cumsum(u):
    if len(u) == 1:
        return u
    elif len(u) == 2:
        return [u[0], u[0]+u[1]]
    else:
        return [u[0]] + cumsum([u[0]+u[1]] + u[2:])

def cumsum_loop(u):
    v = [u[0]]
    for i in u[1:]:
        v.append(v[-1] + i)
    return v

def parse(fmt, s):
    """ Based on an ASCII format string, parse the information in *s*. """
    nch = reclen(fmt)
    reps = nreps(fmt)
    t = dtype(fmt)
    pos = [0] + cumsum(nch*reps)
    vals = [t_(s[a:b]) for a,b,t_ in zip(pos[:-1], pos[1:],
                                             itertools.cycle(t))]
    return vals if len(vals) > 1 else vals[0]


