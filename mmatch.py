#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function


# =============================================================================
# DOCS
# =============================================================================

__doc__ = """foo"""

__version__ = "0.0.1"


# =============================================================================
# IMPORTS
# =============================================================================

import sys
import os
import shutil
import argparse
import logging
import multiprocessing as mp
import csv
import glob
import tempfile
import atexit
import uuid
import copy

import sh

import numpy as np


# =============================================================================
# LOG
# =============================================================================

logger = logging.getLogger("mmatch")
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.WARNING)


# =============================================================================
# CONSTANTS
# =============================================================================

DOC = __doc__

VERSION = __version__

DEFAULT_RADIUS = 3 * 9.2592592592592588e-5

try:
    DEFAULT_VVV_FLX2MAG = sh.Command("./vvv_flx2mag")
except sh.CommandNotFound:
    DEFAULT_VVV_FLX2MAG = None
    logger.warning("Default 'vvv_flx2mag' not found")

CPU_COUNT = mp.cpu_count()

SOURCE_DTYPE = {
    "names": ['ra_h', 'ra_m', 'ra_s', 'dec_d', 'dec_m', 'dec_s'],
    "formats": [int, int, float, int, int, float]
}

PAWPRINT_DTYPE = {
    "names": [
        'ra_h', 'ra_m', 'ra_s', 'dec_d', 'dec_m', 'dec_s', 'x', 'y',
        'mag1', 'mag_err1', 'mag2', 'mag_err2',
        'mag3', 'mag_err3', 'mag4', 'mag_err4',
        'mag5', 'mag_err5', 'mag6', 'mag_err6', 'mag7', 'mag_err7',
        'chip_nro', 'stel_cls', 'elip', 'pos_ang', 'confidence',
    ],
    "formats": [
        int, int, float, int, int, float, float, float,
        float, float, float, float,
        float, float, float, float,
        float, float, float, float, float, float,
        int, int, float, float, float
    ]
}



# =============================================================================
# CLASSES
# =============================================================================

class Matcher(mp.Process):

    def __init__(self, idx, ras, decs, pawprints, vvv_flx2mag,
                 work_directory, radius):
        super(Matcher, self).__init__()
        self.idx = idx
        self.ras = ras
        self.decs = decs
        self.pawprints = pawprints
        self.vvv_flx2mag = vvv_flx2mag
        self.radius = radius
        self.uuid = str(uuid.uuid1())
        self.work_directory = work_directory
        self.ascii_directory = os.path.join(work_directory, "ascii")
        self.array_directory = os.path.join(work_directory, "arrays")
        self.setup_dirs()

    def setup_dirs(self):
        paths = (self.ascii_directory,
                 self.array_directory)
        for path in paths:
            if not os.path.exists(path):
                os.makedirs(path)

    def load_pawprint(self, pawprint):
        basename = os.path.splitext(os.path.basename(pawprint))[0]
        arrayname = "{}.npy".format(basename)
        arraypath = os.path.join(self.array_directory, arrayname)
        if not os.path.exists(arraypath):
            asciiname = "{}.txt".format(basename)
            asciipath = os.path.join(self.ascii_directory, asciiname)
            if not os.path.exists(asciipath):
                # create the ascii table
                self.vvv_flx2mag(pawprint, asciipath)

            # read ascii table
            odata = np.genfromtxt(asciipath, PAWPRINT_DTYPE)

            # extract the ra and dec as degrees
            radeg, decdeg = radec_deg(odata)

            # create a new dtype to store the ra and dec as degrees
            dtype = copy.deepcopy(PAWPRINT_DTYPE)
            dtype["names"].insert(0, "dec_deg")
            dtype["names"].insert(0, "ra_deg")
            dtype["formats"].insert(0, float)
            dtype["formats"].insert(0, float)

            # create an empty array and copy the values
            data = np.empty(len(odata), dtype=dtype)
            for name in data.dtype.names:
                if name == "ra_deg":
                    data[name] = radeg
                elif name == "dec_deg":
                    data[name] = decdeg
                else:
                    data[name] = odata[name]

            # store the numpy array to future uses
            np.save(arraypath, data)
            return data

        return np.load(arraypath)

    def run(self):
        logger.info("Starting matcher '{}'...".format(self.idx))

        tile_ra, tile_dec = self.ras, self.decs

        for pawprint in self.pawprints:
            pwp_data = self.load_pawprint(pawprint)

            pwp_ra, pwp_dec = pwp_data["ra_deg"], pwp_data["dec_deg"]

            import ipdb; ipdb.set_trace()

        # match the sources
        # store the value
        logger.info("Matcher '{}[]' DONE!".format(self.idx, self.uuid))


# =============================================================================
# FUNCTIONS
# =============================================================================

def chunk_it(seq, num):
    """Split a sequence in a 'num' sequence of ~same size

    """
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return sorted(out, reverse=True)


def radec_source(string):
    """Convert a RA, DEC string with the format "17:29:21.4 -30:56:02" to
    a numpy record array
    """
    ra, dec = [map(float, e.split(":")) for e in string.split()]
    row = tuple(ra + dec)
    return np.array([row], dtype=SOURCE_DTYPE)


def csv_sources(path):
    """Take the first 2 columns of a CSV and convert it to a
    numpy record array

    """
    ra_decs = []
    with open(path) as fp:
        for row in csv.reader(fp):
            ra_dec_string = " ".join(row[:2])
            ra_decs.append(radec_source(ra_dec_string))
    return ra_decs


def radec_deg(sources):
    """Generate two arrays with RA and DEC as degree

    """
    ra = 15 * (sources['ra_h'] +
               sources['ra_m'] / 60.0 +
               sources['ra_s'] / 3600.0)

    dec = np.sign(sources['dec_d']) * (np.abs(sources['dec_d']) +
                                       sources['dec_m'] / 60.0 +
                                       sources['dec_s'] / 3600.0)

    return ra, dec



# =============================================================================
# MAIN
# =============================================================================

def _main(argv):

    def get_parser():
        parser = argparse.ArgumentParser(description=DOC, version=VERSION)

        src_group = parser.add_mutually_exclusive_group(required=True)
        src_group.add_argument(
            '--sources-file', '-srcf', dest='sources', action='store',
            type=csv_sources, metavar="PATH",
            help='CSV file with at least 2 columns "ra, dec"')
        src_group.add_argument(
            '--sources', '-srcs', dest='sources', action='store',
            type=radec_source, metavar='"RA DEC"', nargs="+",
            help='Quoted ra and dec parameters example: '
                 '"17:29:21.4 -30:56:02"')

        parser.add_argument(
            '--quiet', '-q', dest='loglevel', action='store_const',
            const=logging.WARNING, default=logging.INFO,
            help='set log level to warning')

        parser.add_argument(
            "--pawprints", "-pwps", dest="pawprints", nargs="+",
            action="store", required=True, metavar="PATH", type=str,
            help="Pawprint files or directory")
        parser.add_argument(
            "--radius", "-r", dest="radius", action="store",
            default=DEFAULT_RADIUS, type=float,
            help="Radiuous to make the crossmatch")
        parser.add_argument(
            "--vvv-flx2mag", "-f2m", dest="vvv_flx2mag",
            required=not bool(DEFAULT_VVV_FLX2MAG), metavar="PATH",
            default=DEFAULT_VVV_FLX2MAG, type=sh.Command,
            help="vvv_flx2mag command PATH")
        parser.add_argument(
            "--working-directory", "-wd", dest="work_directory",
            action="store", metavar="PATH",
            help="directory to store the temporary files")
        parser.add_argument(
            '--output', '-o', dest='output', action='store',
            type=argparse.FileType('w'), metavar="PATH",
            help='destination of your pawprint')
        parser.add_argument(
            "--procs", "-p", dest="procs", default=CPU_COUNT,
            metavar="NUMBER", type=int, choices=list(range(CPU_COUNT+1)),
            help="Number of processors to use (0 run in the same process)")
        return parser

    def to_files(paths):
        files = []
        for path in paths:
            abspath = os.path.abspath(path)
            if os.path.isdir(abspath):
                pattern = os.path.join(abspath, "*.fits")
                files.extend([
                    p for p in glob.glob(pattern) if os.path.isfile(p)])
            else:
                files.append(abspath)
        return list(set(files))

    def work_dir(path):
        if path:
            return path
        path = tempfile.mkdtemp(suffix="_mmatch")
        atexit.register(shutil.rmtree, path=path)
        return path

    # parse the arguments
    parser = get_parser()
    args = parser.parse_args(argv)

    # extract and post-process the input data
    src_ras, src_decs = radec_deg(np.concatenate(args.sources))
    loglevel = args.loglevel
    pawprints = to_files(args.pawprints)
    radius = args.radius
    vvv_flx2mag = args.vvv_flx2mag
    work_directory = work_dir(args.work_directory)
    output = args.output
    procs = args.procs

    logger.setLevel(loglevel)
    logger.info("Using vvv_flx2mag: {}".format(vvv_flx2mag))

    # at least we need one proc to run
    if procs:
        procs_to_span = procs
        logger.info("Starting {} matchers".format(procs))
    else:
        logger.info("Synchronous match")
        procs_to_span = 1

    # run the process
    running_procs = []
    for idx, chunk in enumerate(chunk_it(pawprints, procs_to_span)):
        matcher = Matcher(
            idx=idx, ras=src_ras, decs=src_decs,
            pawprints=chunk, vvv_flx2mag=vvv_flx2mag,
            work_directory=work_directory, radius=radius)
        if procs:
            matcher.start()
            running_procs.append(matcher)
        else:
            matcher.run()
    for matcher in running_procs:
        matcher.join()

    # combine the outputs


if __name__ == "__main__":
    _main(sys.argv[1:])
