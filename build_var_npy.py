#!/usr/bin/env python3
"""Imputation challenge scoring script

Author:
    Jin Lee (leepc12@gmail.com)
"""

import numpy
import pyBigWig
from score_metrics import normalize_dict
from bw_to_npy import write_dict_to_npy, load_npy
from logger import log


def build_var_dict(npys, chroms):
    y_all = {}
    for c in chroms:
        y_all[c] = []

    for f in npys:
        y_dict = load_npy(f)
        y_dict_norm = normalize_dict(y_dict, chroms)

        for c in chroms:
            y_all[c].append(y_dict_norm[c])

    var = {}
    for c in chroms:
        var[c] = numpy.std(numpy.array(y_all[c]), axis=0) ** 2

    return var


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(
        description='ENCODE Imputation Challenge variance .npy builder')
    parser.add_argument('npy', nargs='+',
                        help='Binned truth .npy file')
    parser.add_argument('--out-npy-prefix', required=True,
                         help='Output prefix for .npy or .npz')
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
    args = parser.parse_args()

    # some submission files have whitespace in path...
    for i, f in enumerate(args.npy):
        args.npy[i] = f.strip("'")
    if args.chrom == ['all']:
        args.chrom = ['chr' + str(i) for i in range(1, 23)] + ['chrX']

    return args


def main():
    args = parse_arguments()

    var = build_var_dict(args.npy, args.chrom)
    write_dict_to_npy(var, args.out_npy_prefix)

    log.info('All done')


if __name__ == '__main__':
    main()
