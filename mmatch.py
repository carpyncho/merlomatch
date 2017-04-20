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

import sh

import numpy as np


# =============================================================================
# CONSTANTS AND LOG
# =============================================================================

DOC = __doc__

DEFAULT_RADIUS = 3 * 9.2592592592592588e-5

DEFAULT_WORK_DIRECTORY = "_temp"

CPU_COUNT = mp.cpu_count()

DTYPE = {
    "names": ['ra_h', 'ra_m', 'ra_s', 'dec_d', 'dec_m', 'dec_s'],
    "formats": [int, int, float, int, int, float]
}


# =============================================================================
# CLASSES
# =============================================================================

class Matcher(mp.Process):

    def __init__(self):
        pass

    def ra_to_degree(self, arr):
        return 15 * (
            arr['ra_h'] +
            arr['ra_m'] / 60.0 +
            arr['ra_s'] / 3600.0)

    def dec_to_degree(self, arr):
        return np.sign(arr['dec_d']) * (
            np.abs(arr['dec_d']) +
            arr['dec_m'] / 60.0 +
            arr['dec_s'] / 3600.0)


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


def ra_dec(string):
    """Convert a RA, DEC string with the format "17:29:21.4 -30:56:02" to
    a numpy record array
    """
    ra, dec = [map(float, e.split(":")) for e in string.split()]
    row = tuple(ra + dec)
    return np.array([row], dtype=DTYPE)


def csv_ra_decs(path):
    """Take the first 2 columns of a CSV and convert it to a
    numpy record array

    """
    ra_decs = []
    with open(path) as fp:
        for row in csv.reader(fp):
            ra_dec_string = " ".join(row[:2])
            ra_decs.append(ra_dec(ra_dec_string))
    return ra_decs


# =============================================================================
# MAIN
# =============================================================================

def _main(argv):

    def get_parser():
        parser = argparse.ArgumentParser(description=DOC)

        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument(
            '--sources-file', '-srcf', dest='sources', action='store',
            type=csv_ra_decs, metavar="PATH",
            help='CSV file with at least 2 columns "ra, dec"')
        group.add_argument(
            '--sources', '-srcs', dest='sources', action='store',
            type=ra_dec, metavar='"RA DEC"', nargs="+",
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
            "-w", "--working-directory", dest="work_directory",
            default=DEFAULT_WORK_DIRECTORY, action="store", metavar="PATH",
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

    # parse the arguments
    parser = get_parser()
    args = parser.parse_args(argv)

    # extract and post-process the input data
    sources = np.concatenate(args.sources)
    pawprints = to_files(args.pawprints)
    radius = args.radius
    work_directory = args.work_directory
    output = args.output
    procs = args.procs

    # at least we need one proc to run
    procs_to_span = procs if procs else 1

    # run the process
    running_procs = []
    for idx, chunk in enumerate(chunk_it(pawprints, procs_to_span)):
        matcher = Match(
            idx=idx, sources=sources, pawprint=chunk,
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
