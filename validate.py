# Imputed track evaluations
# Author: Jacob Schreiber, Jin Lee
# Contact: jmschreiber91@gmail.com, leepc12@gmail.com

import sys
import argparse
import pyBigWig
import json
import logging

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)

# cat hg38.chrom.sizes | grep -P "chr[\dX]" | grep -v _
CHRSZ = {
    'chr1': 248956422,
    'chr2': 242193529,
    'chr3': 198295559,
    'chr4': 190214555,
    'chr5': 181538259,
    'chr6': 170805979,
    'chr7': 159345973,
    'chr8': 145138636,
    'chr9': 138394717,
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
    'chr20': 64444167,
    'chr21': 46709983,
    'chr22': 50818468,
    'chrX': 156040895
}

WINDOW_SIZE = 25


def validate(bw):
    # print chrsz
    log.info('Validating chromosome sizes...')
    log.info('==== Template ====')
    print(json.dumps(CHRSZ, indent=4))
    log.info('==== Your submission ====')
    print(json.dumps(bw.chroms(), indent=4))
    # check number of chrs
    if len(bw.chroms()) != len(CHRSZ):
        raise ValueError('Number of chrs {}. It must be {}'.format(
            len(bw.chroms()), len(CHRSZ)))
    # check each chrsz
    for k, v in bw.chroms().items():
        if k not in CHRSZ:
            raise ValueError('Invalid chr {}'.format(k))
        if v != CHRSZ[k]:
            raise ValueError('Invalid size {} for chr {}'.format(v, k))

    # validate window_size
    for c, s in CHRSZ.items():
        log.info('Validating chromosome: {}...'.format(c))
        for start, end, v in bw.intervals(c):
            if end == s:
                continue
            if end-start != WINDOW_SIZE:
                raise ValueError(
                    'Invalid window size for chr {}. '
                    'start: {}, end: {}, value: {}'.format(
                        c, start, end, v))
    log.info('Validation done.')
    return 0


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

    validate(bw)

    return 0


if __name__ == '__main__':
    main()
