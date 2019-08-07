#!/usr/bin/env python3
"""Imputation challenge submission bigwig validator

Author:
    Jin Lee (leepc12@gmail.com)
"""

import pyBigWig
import json
import traceback
import math
from logger import log


# cat hg38.chrom.sizes | grep -P "chr[\dX]" | grep -v _
CHRSZ = {
    'chr1': 248956422,
    'chr10': 133797422,
    'chr11': 135086622,
    'chr12': 133275309,
    'chr13': 114364328,
    'chr14': 107043718,
    'chr15': 101991189,
    'chr16': 90338345,
    'chr17': 83257441,
    'chr18': 80373285,
    'chr19': 58617616,
    'chr2': 242193529,
    'chr20': 64444167,
    'chr21': 46709983,
    'chr22': 50818468,
    'chr3': 198295559,
    'chr4': 190214555,
    'chr5': 181538259,
    'chr6': 170805979,
    'chr7': 159345973,
    'chr8': 145138636,
    'chr9': 138394717,
    'chrX': 156040895
}

def validate(bw_file, window_size=25):
    log.info('Opening bigwig file...')
    bw = pyBigWig.open(bw_file.strip("'"))
    all_msg = ''
    
    try:
        valid = True
        # print chrsz
        # log.info('Validating chromosome sizes...')
        # log.info('==== Template ====')
        # print(json.dumps({k: CHRSZ[k] for k in sorted(CHRSZ)}, indent=4))
        # log.info('==== Your submission ====')
        # print(json.dumps({k: bw.chroms()[k] for k in sorted(bw.chroms())}, indent=4))
        # check number of chrs
        if len(bw.chroms()) != len(CHRSZ):
            msg = 'Invalid number of chromosome {}. It should match with {}'.format(
                len(bw.chroms()), len(CHRSZ))
            print(msg)
            all_msg += msg + '\n'
            valid = False
        # check each chrsz
        for k, v in bw.chroms().items():
            if k not in CHRSZ:
                msg = 'Invalid chromosome {}'.format(k)
                print(msg)
                all_msg += msg + '\n'
                valid = False
                continue
            if v != CHRSZ[k]:
                msg = 'Invalid size {} for chromosome {}'.format(v, k)
                print(msg)
                all_msg += msg + '\n'
                valid = False
                continue

        # validate window_size
        for c, s in CHRSZ.items():
            log.info('Validating chromosome {}...'.format(c))
            if bw.intervals(c) is None:
                msg = 'No intervals found for chromosome {}. '.format(c)
                print(msg)
                all_msg += msg + '\n'
                valid = False
                continue
            for start, end, v in bw.intervals(c):
                if end == s:
                    continue
                if end-start != window_size:
                    msg = 'Invalid window size for chromosome {}. '\
                    'start: {}, end: {}, value: {}'.format(
                            c, start, end, v)
                    print(msg)
                    all_msg += msg + '\n'
                    valid = False
                if end > s:
                    msg = 'Invalid end interval in chromosome {}. '\
                    'End must be equal to or smaller than chrom size. '
                    'start: {}, end: {}, value: {}, chrsz: {}'.format(
                            c, start, end, v, s)
                    print(msg)
                    all_msg += msg + '\n'
                    valid = False
                if math.isnan(v) or v == float('inf') or v == float('-inf'):
                    msg = 'Found NaN or Inf. '\
                    'start: {}, end: {}, value: {}, chrsz: {}'.format(
                            c, start, end, v, s)
                    print(msg)
                    all_msg += msg + '\n'
                    valid = False

    except Exception as e:
        st = StringIO()
        traceback.print_exc(file=st)
        msg = st.getvalue()
        print(msg)
        all_msg += msg + '\n'
        valid = False

    if valid:
        log.info('Validation done successfully.')
        return True, all_msg
    else:
        log.info('Validation failed.')
        return False, all_msg

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description='ENCODE Imputation Challenge'
                                          'validation script.')
    parser.add_argument('bw', type=str,
                        help='Bigwig file to be validated.')
    parser.add_argument('--window-size', default=25, type=int,
                         help='Window size for bigwig in bp')
    args = parser.parse_args()
    return args


def main():
    # read params
    args = parse_arguments()

    valid, _ = validate(args.bw, args.window_size)
    if valid:
        return 0
    else:
        return 1

if __name__ == '__main__':
    main()
