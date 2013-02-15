""" This module includes the `GeoJSON` class and convenience functions for
dealing with GeoJSON data. The `GeoJSON` class uses the builtin json module and
enforces GeoJSON standards. """

import json
import traceback

__status__ = 'incomplete'

def geojson2csv(fin, fout):
    """ Convert a JSON file to a CSV. Only properties (attributes) common to
    all features are retained. The arguments `fin` and `fout` may be either
    filenames or file-like objects. """

    if hasattr(fin, 'read'):
        jsdat = json.load(fin)
    else:
        with open(fin, 'r') as f:
            jsdat = json.load(f)

    # Load coordinates and properties
    X = []
    Y = []
    PROP = []
    try:
        for feat in jsdat['features']:
            x, y = feat['geometry']['coordinates']
            X.append(x)
            Y.append(y)
            PROP.append(feat['properties'])
    except KeyError:
        sys.stderr.write('invalid geojson schema')
        traceback.print_exc()
    except:
        traceback.print_exc()

    # Find subset of properties that exist in all features
    if len(PROP) > 1:
        s = set(PROP[0].keys())
        common_props = s.intersection(*[set(p.keys()) for p in PROP[1:]])
    else:
        common_props = PROP[0].keys()


    # Output
    outstr = 'x,y' + reduce(lambda a,b: a+','+b, common_props) + '\n'   # Construct header

    outstr += reduce(lambda a, b: a+'\n'+b,                         # Combine lines
                map(lambda x, y, p: str(x) + ',' + str(y) + ',' +   # Construct lines
                    reduce(lambda a, b: str(a) + ',' + str(b),      # Construct properties
                        [p[cp] for cp in common_props]),            # Filter spurious props
                    X, Y, PROP))                                    # Args for line lambda

    try:
        if hasattr(fout, 'read'):
            f = fout
        else:
            f = open(fout, 'w')
    
        f.write(outstr)

    finally:
        if not hasattr(fout, 'read'):
            f.close()
    
    return


