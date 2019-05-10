# Imputed track evaluations
# Author: Jacob Schreiber, Jin Lee
# Contact: jmschreiber91@gmail.com, leepc12@gmail.com

import sys
import argparse
import pyBigWig
import json
import logging
import traceback

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)

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

WINDOW_SIZE = 25


def validate(bw):
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
                if end-start != WINDOW_SIZE:
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
    parser = argparse.ArgumentParser(prog='ENCODE Imputation Challenge'
                                          'validation script.')
    parser.add_argument('bw', type=str,
                        help='Bigwig file to be validated.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING',
                                 'CRITICAL', 'ERROR', 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def main():
    # read params
    args = parse_arguments()

    log.info('Opening bigwig file...')
    bw = pyBigWig.open(args.bw)

    return validate(bw)


if __name__ == '__main__':
    main()
