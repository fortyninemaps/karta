""" Low-level functions for reading ESRI ASCII grids """

import numpy as np

def aairead(fnm):
    """ Read an existing ASCII grid file and return a Numpy array and a
    dictionary of header information. """
    try:
        h = []
        with open(fnm, 'r') as f:
            # Count the number of header entries
            cnt = 0
            while True:
                l = f.readline()
                if l.split(None, 1)[0].lower() in ['nrows', 'ncols',
                        'yllcenter', 'xllcenter', 'yllcorner', 'xllcorner',
                        'cellsize', 'nodata_value']:
                    cnt += 1
                else:
                    break

            # Read the header, then the array data
            f.seek(0)
            for _ in range(cnt):
                h.append(f.readline())
            data = f.readlines()

    except IOError:
        raise AAIIOError('error while trying to open {0}'.format(fnm))

    # Store data internally
    hs = [rec.split(None, 1) for rec in h]
    hdr = dict([(rec[0].lower(), float(rec[1])) for rec in hs])
    hdr['ncols'] = int(hdr['ncols'])
    hdr['nrows'] = int(hdr['nrows'])
    if 'yllcenter' not in hdr.keys():
        hdr['yllcenter'] = None

    if 'xllcenter' not in hdr.keys():
        hdr['xllcenter'] = None

    if 'yllcorner' not in hdr.keys():
        hdr['yllcorner'] = None

    if 'xllcorner' not in hdr.keys():
        hdr['xllcorner'] = None

    if 'nodata_value' not in hdr.keys():
        hdr['nodata_value'] = -9999

    check_header(hdr)

    f = lambda l: [float(i) for i in l.split()]
    data_f = list(map(f, data))
    data_a = np.array(data_f)
    data_a[data_a==hdr['nodata_value']] = np.nan

    return data_a, hdr


def check_header(hdr):
    """ Make sure that all required header records are present, as well as both
    centered and corner spatial references. """
    for field in ('ncols', 'nrows', 'cellsize'):
        if field not in hdr:
            raise AAIError("{0} not defined".format(field.upper()))

    d = hdr['cellsize']

    if hdr.get('yllcenter') is None:
        if hdr.get('yllcorner') is None:
            raise AAIIOError('YLL reference not defined')
        else:
            hdr['yllcenter'] = hdr['yllcorner'] + d / 2.0
    else:
        hdr['yllcorner'] = hdr['yllcenter'] - d / 2.0

    if hdr.get('xllcenter') is None:
        if hdr.get('xllcorner') is None:
            raise AAIIOError('XLL reference not defined')
        else:
            hdr['xllcenter'] = hdr['xllcorner'] + d / 2.0
    else:
        hdr['xllcorner'] = hdr['xllcenter'] - d / 2.0

    return hdr


class AAIError(Exception):
    def __init__(self, message=''):
        self.message = message
    def __str__(self):
        return self.message


class AAIIOError(Exception):
    """ Exceptions related to ESRI ASCII grid driver. """
    def __init__(self, value, detail=None):
        self.value = value
        self.detail = detail
    def __str__(self):
        if self.detail is not None:
            s = self.value + ': ' + self.detail
        else:
            s = self.value
        return repr(s)
