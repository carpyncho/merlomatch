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
import uuid

import sh

import numpy as np

from astropysics import coords


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

MATCH_FIELDS = (
    ["src.idx"]  + ["src.{}".format(n) for n in SOURCE_DTYPE["names"]] +
    ["pwp.path"] + ["pwp.{}".format(n) for n in PAWPRINT_DTYPE["names"]])


# =============================================================================
# CLASSES
# =============================================================================

class Matcher(mp.Process):

    def __init__(self, idx, sources, pawprints, vvv_flx2mag,
                 work_directory, radius):
        super(Matcher, self).__init__()
        self.idx = idx
        self.sources = sources
        self.pawprints = pawprints
        self.vvv_flx2mag = vvv_flx2mag
        self.radius = radius
        self.uuid = str(uuid.uuid1())
        self.work_directory = work_directory
        self.ascii_directory = os.path.join(work_directory, "ascii")
        self.array_directory = os.path.join(work_directory, "arrays")
        self.temp_directory = os.path.join(work_directory, "temp")
        self.tempfile = "{}.csv".format(self.uuid)
        self.tempfile_path = os.path.join(self.temp_directory, self.tempfile)
        self.setup_dirs()

    def setup_dirs(self):
        paths = (self.ascii_directory,
                 self.array_directory,
                 self.temp_directory)
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
            data = radec_deg(odata, PAWPRINT_DTYPE)

            # store the numpy array to future uses
            np.save(arraypath, data)
            return data

        return np.load(arraypath)

    def match(self, tile_ra, tile_dec, pwp_ra, pwp_dec):
        nearestind_pwp, distance_pwp, match_pwp = coords.match_coords(
            tile_ra, tile_dec, pwp_ra, pwp_dec,
            eps=self.radius, mode="nearest")
        nearestind_ms, distance_ms, match_ms = coords.match_coords(
            pwp_ra, pwp_dec, tile_ra, tile_dec,
            eps=self.radius, mode="nearest")
        for idx_pwp, idx_ms in enumerate(nearestind_ms):
            if match_ms[idx_pwp] and \
               nearestind_pwp[idx_ms] == idx_pwp \
               and match_pwp[idx_ms]:
                    yield idx_ms, idx_pwp

    def get_output(self):
        with open(self.tempfile_path) as fp:
            return fp.read()

    def run(self):
        logger.info("Starting matcher '{}'".format(self.idx))
        ras, decs = self.sources["ra_deg"], self.sources["dec_deg"]
        with open(self.tempfile_path, "w") as fp:
            writer = csv.DictWriter(fp, MATCH_FIELDS)
            for pawprint in self.pawprints:
                logger.info("Matcher {}: '{}'...".format(self.idx, pawprint))

                pwp_data = self.load_pawprint(pawprint)
                pwp_ra, pwp_dec = pwp_data["ra_deg"], pwp_data["dec_deg"]

                matchs = self.match(ras, decs, pwp_ra, pwp_dec)

                for src_idx, pwp_idx in matchs:
                    source, pwp = self.sources[src_idx], pwp_data[pwp_idx]
                    row = {"src.idx": src_idx, "pwp.path": pawprint}
                    for field in MATCH_FIELDS:
                        if field in row:
                            continue
                        elif field.startswith("pwp."):
                            row[field] = pwp[field.split(".")[-1]]
                        else:
                            row[field] = source[field.split(".")[-1]]
                    writer.writerow(row)
        logger.info("Matcher '{}' DONE!".format(self.idx))


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


def add_deg_columns(odata, dtypes, radeg, decdeg):
    """Add ra_deg and dec_deg columns to existing recarray

    """
    # create a new dtype to store the ra and dec as degrees
    dtype = copy.deepcopy(dtypes)
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
    return data


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


def radec_deg(sources, dtypes):
    """Generate two arrays with RA and DEC as degree

    """
    ra = 15 * (sources['ra_h'] +
               sources['ra_m'] / 60.0 +
               sources['ra_s'] / 3600.0)

    dec = np.sign(sources['dec_d']) * (np.abs(sources['dec_d']) +
                                       sources['dec_m'] / 60.0 +
                                       sources['dec_s'] / 3600.0)

    return add_deg_columns(sources, dtypes, ra, dec)



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
    sources = radec_deg(np.concatenate(args.sources), SOURCE_DTYPE)
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
            idx=idx, sources=sources,
            pawprints=chunk, vvv_flx2mag=vvv_flx2mag,
            work_directory=work_directory, radius=radius)
        if procs:
            matcher.start()
        else:
            matcher.run()
        running_procs.append(matcher)

    logger.info("Mergin outputs...")
    output.write(",".join(MATCH_FIELDS) + "\n")
    for matcher in running_procs:
        if matcher.is_alive():
            matcher.join()
        output.write(matcher.get_output())


if __name__ == "__main__":
    _main(sys.argv[1:])
