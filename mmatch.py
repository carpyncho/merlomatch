#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, unicode_literals


# =============================================================================
# DOCS
# =============================================================================

__doc__ = """foo"""

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

import sh

import numpy as np


# =============================================================================
# LOG
# =============================================================================

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("mmatch")


# =============================================================================
# CONSTANTS
# =============================================================================

DOC = __doc__

DEFAULT_RADIUS = 3 * 9.2592592592592588e-5

try:
    DEFAULT_FITSIO_CAT_LIST = sh.Command("./fitsio_cat_list")
except sh.CommandNotFound:
    DEFAULT_FITSIO_CAT_LIST = None
    logger.warning("Default 'fitsio_cat_list' not found")

CPU_COUNT = mp.cpu_count()

DTYPE = {
    "names": ['ra_h', 'ra_m', 'ra_s', 'dec_d', 'dec_m', 'dec_s'],
    "formats": [int, int, float, int, int, float]
}

# =============================================================================
# CLASSES
# =============================================================================

class Matcher(mp.Process):

    def __init__(self, idx, ras, decs, pawprints, fitsio_cat_list,
                 work_directory, radius):
        super(Matcher, self).__init__()
        self.idx = idx
        self.ras = ras
        self.decs = decs
        self.pawprints = pawprints
        self.fitsio_cat_list = fitsio_cat_list
        self.radius = radius
        self.uuid = str(uuid.uuid1())
        self.work_directory = work_directory
        self.ascii_directory = os.path.join(work_directory, "ascii")
        self.array_directory = os.path.join(work_directory, "arrays")
        self.temp_directory_directory = os.path.join(
            work_directory, "temp", self.uuid)
        self.setup_dirs()

    def setup_dirs(self):
        paths = (self.ascii_directory,
                 self.array_directory, self.temp_directory_directory)
        for path in paths:
            if not os.path.exists(path):
                os.makedirs(path)

    def load_pawprint(self, pawprint):
        basename = os.path.splitext(os.path.basename(pawprint))[0]
        arrayname = "{}.npy".format(basename)
        arraypath = os.path.join(self.array_directory, arrayname)
        if not os.path.exists(arraypath):
            asciiname = "{}.txt".format(basename)
            asciipath = os.path.join(self.ascii_directory, arrayname)
            if not os.path.exists(asciipath):
                # convertir a ascii
            # leer ascii
            # crear el numpy

        return np.load(arraypath)



    def run(self):
        logger.info("Starting matcher '{}[{}]'...".format(self.idx, self.uuid))

        # setup the directory



        for pawprint in self.pawprints:
            pwp_ras, pwp_decs = self.load_pawprint(pawprint)
            # convert to ascii


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
    return np.array([row], dtype=DTYPE)


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
        parser = argparse.ArgumentParser(description=DOC)

        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument(
            '--sources-file', '-srcf', dest='sources', action='store',
            type=csv_sources, metavar="PATH",
            help='CSV file with at least 2 columns "ra, dec"')
        group.add_argument(
            '--sources', '-srcs', dest='sources', action='store',
            type=radec_source, metavar='"RA DEC"', nargs="+",
            help='Quoted ra and dec parameters example: '
                 '"17:29:21.4 -30:56:02"')
        parser.add_argument(
            "--pawprints", "-pwps", dest="pawprints", nargs="+",
            action="store", required=True, metavar="PATH", type=str,
            help="Pawprint files or directory")
        parser.add_argument(
            "--radius", "-r", dest="radius", action="store",
            default=DEFAULT_RADIUS, type=float,
            help="Radiuous to make the crossmatch")
        parser.add_argument(
            "--fitsio-cat-list", "-fcl", dest="fitsio_cat_list",
            required=not bool(DEFAULT_FITSIO_CAT_LIST), metavar="PATH",
            default=DEFAULT_FITSIO_CAT_LIST, type=sh.Command,
            help="fitsio_cat_list command PATH")
        parser.add_argument(
            "-w", "--working-directory", dest="work_directory",
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
    pawprints = to_files(args.pawprints)
    radius = args.radius
    fitsio_cat_list = args.fitsio_cat_list
    work_directory = work_dir(args.work_directory)
    output = args.output
    procs = args.procs

    logger.info("Using fitsio_cat_list: {}".format(fitsio_cat_list))

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
            pawprints=chunk, fitsio_cat_list=fitsio_cat_list,
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
