# Bootrapping code
# Author: Jin Lee
# Contact: leepc12@gmail.com

import sys
import numpy
import pyBigWig
import argparse
import logging
from sklearn.model_selection import StratifiedKFold


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE Imputation Challenge bootstrapped label generator')
    parser.add_argument('submission_template_bw',
                        help='Submission template bigwig file.')
    parser.add_argument('--fold', default=10, type=int,
                        help='Fold for bootstrapping '
                             '(n_splits for StratifiedKFold)')
    parser.add_argument('--window-size', default=25,
                        help='Window size for bigwig in bp.')
    parser.add_argument('--random-seed', default=0,
                        help='Random seed (random_state for StratifiedKFold)')
    parser.add_argument('--out-npy', default='bootstrapped_label',
                        help='Write bootstrapped label to .npy file.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING',
                                 'CRITICAL', 'ERROR', 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


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

    # result = [dict()]*args.fold  # DANGER! IT DOES NOT WORK
    result = []
    for k in range(args.fold):
        result.append(dict())

    sfk = StratifiedKFold(n_splits=args.fold,
                          shuffle=True,
                          random_state=args.random_seed)
    for c, chrom_len in bw.chroms().items():
        # number of bins per chromosome
        print(c)
        n_bin = (chrom_len-1)//args.window_size+1
        for bootstrap_id, (label, _) in enumerate(sfk.split(
                                                numpy.zeros(n_bin),
                                                numpy.zeros(n_bin))):
            result[bootstrap_id][c] = label
            # print(bootstrap_id, result[bootstrap_id][c][0:20])

        for k in range(args.fold):
            print('k={}: {}'.format(k, result[k][c][0:20]))

    numpy.save(args.out_npy, result)


if __name__ == '__main__':
    main()
