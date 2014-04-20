""" Functions for reading the USGS DEM format. """

import math
import itertools
import re

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

BLOCKB = [("rc_id",0,12,"2I6"),
          ("mn",12,24,"2I6"),
          ("xy",24,72,"2D24.15"),
          ("zdatum",72,96,"D24.15"),
          ("zminmax",96,144,"2D24.15")]

BLOCKC = [("filestatcode",0,6,"I6"),
          ("file_rmse",6,24,"3I6"),
          ("file_nsmp",24,30,"I6"),
          ("demstatcode",30,36,"I6"),
          ("dem_rmse",36,54,"3I6"),
          ("dem_nsmp",54,60,"I6")]

def nreps(fmt):
    """ Return the number of characters in a single record of
    a potentially multiple-record string. """
    try:
        reps = int(fmt[0])
        return reps * nreps(fmt[1:].lstrip("(").rstrip(")"))

    except ValueError:
        return 1
        #return fmt.count(",") + 1

def nreps_re(fmt):
    M = re.match(r"\A\d+", fmt)
    if M:
        return int(fmt[:M.end()]) * nreps_re(fmt[M.end():].lstrip("(").rstrip(")"))
    else:
        return 1

def reclen(fmt):
    """ Return the number of characters in a single record of
    a potentially multiple-record string. """
    try:
        reps = int(fmt[0])
        return reclen(fmt[1:].lstrip("(").rstrip(")"))

    except ValueError:
        if "," in fmt:
            return list(map(lambda a: reclen(a)[0], fmt.split(",")))

        if fmt[0] in ("G", "F", "E", "D"):
            nch = int(fmt[1:].split(".")[0])
        elif fmt[0] in ("A", "I"):
            nch = int(fmt[1:])
        else:
            raise ValueError("Unrecognized format code {0}"
                             .format(fmt))
        return [nch]

def dtype(fmt):
    """ Return a constructor for the datatype encoded in *fmt*. """
    try:
        reps = int(fmt[0])
        return dtype(fmt[1:].lstrip("(").rstrip(")"))

    except ValueError:
        if "," in fmt:
            return list(map(lambda a: dtype(a)[0], fmt.split(",")))
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
    reps = nreps_re(fmt)
    t = dtype(fmt)
    pos = [0] + cumsum(nch*reps)
    vals = [t_(s[a:b]) for a,b,t_ in zip(pos[:-1], pos[1:],
                                             itertools.cycle(t))]
    return vals if len(vals) > 1 else vals[0]

def read_block(blk, KEY):
    """ Read a block from a USGS .dem file given *KEY*, which is a list of
    tuples in (name::str, pos::int, end::int, format::str) form. """
    if len(blk) > 1024:
        blk = blk[:1024]
    blkdict = {}
    for field,pos0,pos1,fmt in KEY[:20]:
        blkdict[field] = parse(fmt, blk[pos0:pos1])
    return blkdict

def demread(fnm):
    """ Read a USGS (or CDED) .dem file and return a dictionary for block A, a list of raw data values, and a dictionary for block B. """
    with open(fnm, "r") as f:
        data = f.read()

    blocka = read_block(data[:1024], BLOCKA)

    dem = []
    profnz = parse("2I6", data[1036:1048])[0]
    profnch = profnz * 6
    for profnum in range(blocka["size"][1]):
        i = (int(math.ceil((profnch-146) / 1024.)) + 1) * profnum + 1
        bhdr = read_block(data[i*1024:(i + 1) * 1024], BLOCKB)
        #profdem = read_block(data[1168:2048], [("data",0,876,"146I6")])["data"]
        profdem = read_block(data[1024*i:1024*(i+1)],
                             [("data",144,1020,"146I6")])["data"]
        blknz = [170 for _ in range((profnz-146)//170)] + [(profnz-146) % 170]
        for j, nz in enumerate(blknz):
            fmt = str(nz) + "I6"
            profdem.extend(read_block(data[1024*(j+i+1) : 1024*(j+i+2)],
                                      [("data",0,1020,fmt)])["data"])
        dem.extend(profdem)

    blockc = read_block(data[-1024:], BLOCKC)
    return blocka, dem, blockc

