#!/usr/bin/env python3
"""Imputation challenge submission bigwig validator

Author:
    Jin Lee (leepc12@gmail.com)
"""

import pyBigWig
import json
import traceback
from logger import log
from challenge_metadata import CHRSZ


def validate(bw_file, window_size=25):
    log.info('Opening bigwig file...')
    bw = pyBigWig.open(bw_file.strip("'"))

    try:
        valid = True
        # print chrsz
        log.info('Validating chromosome sizes...')
        log.info('==== Template ====')
        print(json.dumps({k: CHRSZ[k] for k in sorted(CHRSZ)}, indent=4))
        log.info('==== Your submission ====')
        print(json.dumps({k: bw.chroms()[k] for k in sorted(bw.chroms())}, indent=4))
        # check number of chrs
        if len(bw.chroms()) != len(CHRSZ):
            print('Invalid number of chromosome {}. It should match with {}'.format(
                len(bw.chroms()), len(CHRSZ)))
            valid = False
        # check each chrsz
        for k, v in bw.chroms().items():
            if k not in CHRSZ:
                print('Invalid chromosome {}'.format(k))
                valid = False
                continue
            if v != CHRSZ[k]:
                print('Invalid size {} for chromosome {}'.format(v, k))
                valid = False
                continue

        # validate window_size
        for c, s in CHRSZ.items():
            log.info('Validating chromosome {}...'.format(c))
            if bw.intervals(c) is None:
                print('No intervals found for chromosome {}. '.format(c))
                valid = False
                continue
            for start, end, v in bw.intervals(c):
                if end == s:
                    continue
                if end-start != window_size:
                    print('Invalid window size for chromosome {}. '
                          'start: {}, end: {}, value: {}'.format(
                            c, start, end, v))
                    valid = False
    except Exception as e:
        traceback.print_exc()
        valid = False

    if valid:
        log.info('Validation done successfully.')
        return 0
    else:
        log.info('Validation failed.')
        return 1

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

    return validate(args.bw, args.window_size)


if __name__ == '__main__':
    main()
