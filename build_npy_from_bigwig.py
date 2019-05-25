#!/usr/bin/env python3
"""Imputation challenge scoring script
Author:
    Jin Lee (leepc12@gmail.com) and Jacob Schreiber (jmschreiber91@gmail.com)
"""

import sys
import os
import argparse
import numpy
import pyBigWig
import logging
from score import read_annotation_bed

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)


def parse_arguments():
    py_path = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser(
        prog='ENCODE Imputation Challenge scoring script')
    parser.add_argument('bw',
                        help='Bigwig file or .npy file (for blacklist filtering)')
    parser.add_argument('--out-npy-prefix',
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
    parser.add_argument('--blacklist-file',
                         default=os.path.join(
                            py_path,
                            'annot/hg38/hg38.blacklist.bed.gz'),
                         help='Blacklist BED file. Bootstrap label will be '
                              'generated after removing overlapping regions '
                              'defined in this file.')
    p_score.add_argument('--window-size', default=25, type=int,
                         help='Window size for bigwig in bp')
    p_score.add_argument('--validated', action='store_true',
                         help='For validated submissions '
                              'with fixed interval length of 25 and valid '
                              'chromosome lengths. It will skip interpolation')
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


def bw_to_dict(bw, chrs, window_size=25, validated=False):
    """
    Returns:
        { chr: [] } where [] is a numpy 1-dim array
    """
    result = {}
    for c in chrs:
        log_msg = 'Reading chromosome {} from bigwig...'.format(c)
        log.info(log_msg)
        result_per_chr = []
        chrom_len = bw.chroms()[c]

        num_step = (chrom_len-1)//window_size+1

        if validated:
            all_steps = bw.intervals(c)
            assert(num_step==len(all_steps))

            for step in range(num_step):
                start = step*window_size
                end = min((step+1)*window_size, chrom_len)
                result_per_chr.append(all_steps[step][2])
        else:
            # reshape raw vector as (num_step, window_size)
            raw = bw.values(c, 0, chrom_len, numpy=True)
            reshaped = numpy.zeros((num_step*window_size,))
            reshaped[:raw.shape[0]] = raw
            # pyBigWig returns nan for values out of bounds
            # convert nan to zero
            x = numpy.nan_to_num(reshaped)
            y = numpy.reshape(x, (-1, window_size))
            # bin it
            # reduce dimension to (num_step, 0) by averaging
            # all values in a step
            result_per_chr = y.mean(axis=1)

            # special treatment for last step (where the first nan is)
            # above averaging method does not work with the end step
            # bw.intervals(c)[-1] is the last interval in bigwig
            last_step = bw.intervals(c)[-1][1]//window_size
            start = last_step*window_size
            end = min((last_step+1)*window_size, chrom_len)
            stat = bw.stats(c, start, end, exact=True)
            if stat[0] is None:
                result_per_chr[last_step]=0.0
            else:
                result_per_chr[last_step]=stat[0]

            # # too expensive
            # for step in range(num_step):
            #     start = step*window_size
            #     end = min((step+1)*window_size, chrom_len)
            #     stat = bw.stats(c, start, end, exact=True)
            #     if stat[0] is None:
            #         result_per_chr.append(0)
            #     else:
            #         result_per_chr.append(stat[0])

        result[c] = numpy.array(result_per_chr)
    return result


def blacklist_filter(d, blacklist):
    result = {}
    for c in d:
        result_per_chr = d[c]

        # remove bins overlapping blacklisted region
        if blacklist is None:
            bfilt_result_per_chr = result_per_chr
        else:
            bfilt_result_per_chr = []
            for i, val in enumerate(result_per_chr):
                if i in blacklist[c]:
                    continue
                else:
                    bfilt_result_per_chr.append(val)

        result[c] = numpy.array(bfilt_result_per_chr)
    return result


def get_blacklist_bin_ids(blacklist, chroms, window_size=25):
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
        bins[c] |= set(range(start_bin_id, end_bin_id))

    # convert into list and then numpy array
    result = {}
    for c in chroms:
        result[c] = numpy.array(list(bins[c]), dtype=numpy.int64)

    return result


def main():
    # read params
    args = parse_arguments()

    if args.bw.lower().endswith(('bw', 'bigwig')):
        log.info('Opening bigwig file...')
        bw = pyBigWig.open(args.bw)
        y_dict = bw_to_dict(bw, args.chrom, args.window_size,
                            args.validated)
    elif args.bw.lower().endswith(('npy', 'npz')):
        y_dict = numpy.load(args.bw, allow_pickle=True)[()]

    if args.blacklist_file is None:
        bfilt_y_dict = y_dict
    else:
        log.info('Reading from blacklist bed file...')
        blacklist_lines = read_annotation_bed(args.blacklist_file)
        blacklist_bin_ids = get_blacklist_bin_ids(
            blacklist_lines, args.chrom, args.window_size)
        bfilt_y_dict = blacklist_filter(y_dict, blacklist_bin_ids)

    log.info('Writing to npy or npz...')
    if args.out_npy_prefix is None:
        prefix, _ = os.path.splitext(args.bw)
        out_npy_prefix = prefix
    else:
        out_npy_prefix = args.out_npy_prefix

    numpy.save(out_npy_prefix, bfilt_y_dict)

    log.info('All done')


if __name__ == '__main__':
    main()

