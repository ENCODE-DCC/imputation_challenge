#!/usr/bin/env python3
"""Imputation challenge bootstrapped label generation script
Author:
    Jin Lee (leepc12@gmail.com)
"""

import sys
import numpy
import pyBigWig
import argparse
import logging
from sklearn.model_selection import StratifiedKFold
from score import read_annotation_bed


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)


def parse_arguments():
    import os
    py_path = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser(
        prog='ENCODE Imputation Challenge bootstrapped label generator')
    parser.add_argument('submission_template_bw',
                        help='Submission template bigwig file.')
    parser.add_argument('--fold', default=10, type=int,
                        help='Fold for bootstrapping '
                             '(n_splits for StratifiedKFold)')
    parser.add_argument('--window-size', default=25,
                        help='Window size for bigwig in bp.')
    parser.add_argument('--blacklist-file',
                         default=os.path.join(
                            py_path,
                            'annot/hg38/hg38.blacklist.bed.gz'),
                         help='Blacklist BED file. Bootstrap label will be '
                              'generated after removing overlapping regions '
                              'defined in this file.')
    parser.add_argument('--random-seed', default=0,
                        help='Random seed (random_state for StratifiedKFold)')
    parser.add_argument('--out-npy-prefix', default='bootstrapped_label',
                        help='Write bootstrapped label to .npy or .npz file.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING',
                                 'CRITICAL', 'ERROR', 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def get_blacklisted_bins(blacklist, chroms, window_size=25):
    """
    Returns:
        { chrom: [] }: label that overlaps with blackstlisted region
    """
    # make empty sets per chr
    bins = {}
    for c in chroms:
        bins[c] = set()

    for line in blacklist:
        c, start, end = line.split()
        start_bin_id = int(start) // window_size
        end_bin_id = int(end) // window_size + 1
        print(c, start, end, start_bin_id, end_bin_id)

        bins[c] |= set(range(start_bin_id, end_bin_id))

    # convert into list and then numpy array
    result = {}
    for c in chroms:
        result[c] = numpy.array(list(bins[c]), dtype=numpy.int64)

    print(result)
    return result


def main():
    """
    Returns
        [
            { chr: bootstrapped_label }, # 1st fold
            { chr: bootstrapped_label }, # 2nd fold
            ...
        ]
    """
    args = parse_arguments()

    template_bw_file = args.submission_template_bw
    bw = pyBigWig.open(template_bw_file)

    log.info('Reading from blacklist bed file...')
    blacklists = read_annotation_bed(args.blacklist_file)

    blacklisted_bins = get_blacklisted_bins(blacklists, bw.chroms(),
                                            args.window_size)

    # result = [dict()]*args.fold  # DANGER! IT DOES NOT WORK
    result = []
    for k in range(args.fold):
        result.append(dict())

    if args.fold == 1:
        for c, chrom_len in bw.chroms().items():
            # number of bins per chromosome
            print(c)
            blacklist = blacklisted_bins[c]
            n_bin = (chrom_len-1)//args.window_size+1

            org_arr = numpy.array(range(n_bin))
            new_arr = numpy.delete(org_arr, blacklist)

            result[0][c] = new_arr

            for k in range(args.fold):
                print('k={}: {}'.format(k, result[k][c][0:20]))
    else:
        sfk = StratifiedKFold(n_splits=args.fold,
                              shuffle=True,
                              random_state=args.random_seed)
        for c, chrom_len in bw.chroms().items():
            # number of bins per chromosome
            print(c)
            blacklist = blacklisted_bins[c]
            n_bin = (chrom_len-1)//args.window_size+1

            org_arr = numpy.array(range(n_bin))
            new_arr = numpy.delete(org_arr, blacklist)

            n_new = n_bin - len(blacklist)

            for bootstrap_id, (label, _) in enumerate(sfk.split(
                                                    numpy.zeros(n_new),
                                                    numpy.zeros(n_new))):
                result[bootstrap_id][c] = new_arr[label]
                # print(bootstrap_id, result[bootstrap_id][c][0:20])

            for k in range(args.fold):
                print('k={}: {}'.format(k, result[k][c][0:20]))

    numpy.save(args.out_npy_prefix, result)


if __name__ == '__main__':
    main()
