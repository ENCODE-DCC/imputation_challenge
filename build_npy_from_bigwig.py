#!/usr/bin/env python3
"""Imputation challenge scoring script
Author:
    Jin Lee (leepc12@gmail.com) and Jacob Schreiber (jmschreiber91@gmail.com)
"""

import sys
import argparse
import numpy
import pyBigWig
from score import bw_to_dict
import logging

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE Imputation Challenge scoring script')
    parser.add_argument('bw',
                        help='Bigwig file to be converted to .npy')
    parser.add_argument('--out-npy-prefix', required=True,
                         help='Output npy prefix')
    p_score = parser.add_argument_group(
                        title='Scoring parameters')
    p_score.add_argument('--chrom', nargs='+',
                         default=['all'],
                         help='List of chromosomes to be combined to be '
                              'scored. '
                              'Set as "all" (default) to score for all '
                              'chromosomes. '
                              '(e.g. "all" or "chr3 chr21") '
                              'It should be "all" to write scores to DB file')
    p_score.add_argument('--window-size', default=25, type=int,
                         help='Window size for bigwig in bp')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING',
                                 'CRITICAL', 'ERROR', 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    # some submission files have whitespace in path...
    args.bw = args.bw.strip("'")
    if args.chrom == ['all']:
        args.chrom = ['chr' + str(i) for i in range(1, 23)] + ['chrX']

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def main():
    # read params
    args = parse_arguments()

    log.info('Opening bigwig file...')
    bw = pyBigWig.open(args.bw)
    y_dict = bw_to_dict(bw, args.chrom, args.window_size, log)

    log.info('Writing to npy...')
    numpy.save(args.out_npy_prefix, y_dict)

    log.info('All done')


if __name__ == '__main__':
    main()

