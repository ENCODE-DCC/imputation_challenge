#!/usr/bin/env python3
"""Imputation challenge scoring script
Author:
    Jin Lee (leepc12@gmail.com) and Jacob Schreiber (jmschreiber91@gmail.com)
"""

import sys
import argparse
import numpy
import pyBigWig
import os
import gzip
import sqlite3
import time
import gc
from collections import namedtuple
from sklearn.metrics import roc_auc_score
from scipy.stats import norm, spearmanr
import logging

from .build_npy_from_bigwig import read_annotation_bed, bw_to_dict, dict_to_arr

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)


Score = namedtuple(
    'Score',
    ('mse', 'gwcorr', 'gwspear', 'mseprom', 'msegene', 'mseenh',
     'msevar', 'mse1obs', 'mse1imp')
)

ScoreDBRecord = namedtuple(
    'ScoreDBRecord',
    ('submission_id', 'team_id', 'submission_fname', 'cell', 'assay',
     'bootstrap_id') + Score._fields
)

DB_TABLE_SCORE = 'score'
DB_QUERY_INSERT = 'INSERT INTO {table} ({cols}) VALUES ({values});'


def mse(y_true, y_pred):
    return ((y_true - y_pred) ** 2.).mean()


def gwcorr(y_true, y_pred):
    return numpy.corrcoef(y_true, y_pred)[0, 1]


def gwspear(y_true, y_pred):
    return spearmanr(y_true, y_pred)[0]


def mseprom(y_true_dict, y_pred_dict, chroms,
            gene_annotations,
            window_size=25, prom_loc=80):
    """
    Args:
        y_true_dict: truth vector per chromosome
            { chr: y_true } where y_true is a numpy 1-dim array.

        y_pre_dict: predicted vector per chromosome
            { chr: y_pred } where y_pred is a numpy 1-dim array.
    """
    sse, n = 0., 0.

    for chrom in chroms:
        y_true = y_true_dict[chrom]
        y_pred = y_pred_dict[chrom]

        for line in gene_annotations:
            chrom_, start, end, _, _, strand = line.split()
            start = int(start) // window_size
            end = int(end) // window_size + 1

            # if chrom_ in ('chrX', 'chrY', 'chrM'):
            #     continue

            if chrom_ != chrom:
                continue

            if strand == '+':
                sse += ((y_true[start-prom_loc: start] -
                         y_pred[start-prom_loc: start]) ** 2).sum()
                n += y_true[start-prom_loc: start].shape[0]

            else:
                sse += ((y_true[end: end+prom_loc] -
                         y_pred[end: end+prom_loc]) ** 2).sum()
                n += y_true[end: end+prom_loc].shape[0]

    return sse / n


def msegene(y_true_dict, y_pred_dict, chroms,
            gene_annotations,
            window_size=25):
    sse, n = 0., 0.

    for chrom in chroms:
        y_true = y_true_dict[chrom]
        y_pred = y_pred_dict[chrom]

        for line in gene_annotations:
            chrom_, start, end, _, _, strand = line.split()
            start = int(start) // window_size
            end = int(end) // window_size + 1

            # if chrom_ in ('chrX', 'chrY', 'chrM'):
            #     continue

            if chrom_ != chrom:
                continue
            sse += ((y_true[start:end] - y_pred[start:end]) ** 2).sum()
            n += end - start

    return sse / n


def mseenh(y_true_dict, y_pred_dict, chroms,
           enh_annotations,
           window_size=25):
    sse, n = 0., 0.

    for chrom in chroms:
        y_true = y_true_dict[chrom]
        y_pred = y_pred_dict[chrom]

        for line in enh_annotations:
            chrom_, start, end, _, _, _, _, _, _, _, _, _ = line.split()
            start = int(start) // window_size
            end = int(end) // window_size + 1

            if chrom_ != chrom:
                continue
            sse += ((y_true[start:end] - y_pred[start:end]) ** 2).sum()
            n += end - start

    return sse / n


def msevar(y_true, y_pred, y_all=None, var=None):
    """Calculates the MSE weighted by the cross-cell-type variance.

    According to the wiki: Computing this measure involves computing,
    for an assay carried out in cell type x and assay type y, a vector of
    variance values across all assays of type y. The squared error between
    the predicted and true value at each genomic position is multiplied by
    this variance (normalized to sum to 1 across all bins) before being
    averaged across the genome.

    Parameters
    ----------
    y_true: numpy.ndarray, shape=(n_positions,)
        The true signal

    y_pred: numpy.ndarray, shape=(n_positions,)
        The predicted signal

    y_all: numpy.ndarray, shape=(n_celltypes, n_positions)
        The true signal from all the cell types to calculate the variance over.
        mutually exclusive with var

    var: numpy.ndarray, shape=(n_positions,)
        pre-computed var vector
        mutually exclusive with y_all

    Returns
    -------
    mse: float
        The mean-squared error that's weighted at each position by the
        variance.
    """

    if var is None and y_all is None:
        return 0.0
    if var is None:
        var = numpy.std(y_all, axis=0) ** 2

    return ((y_true - y_pred) ** 2).dot(var)/var.sum()


def mse1obs(y_true, y_pred):
    n = int(y_true.shape[0] * 0.01)
    y_true_sorted = numpy.sort(y_true)
    y_true_top1 = y_true_sorted[-n]
    idx = y_true >= y_true_top1

    return mse(y_true[idx], y_pred[idx])


def mse1imp(y_true, y_pred):
    n = int(y_true.shape[0] * 0.01)
    y_pred_sorted = numpy.sort(y_pred)
    y_pred_top1 = y_pred_sorted[-n]
    idx = y_pred >= y_pred_top1

    return mse(y_true[idx], y_pred[idx])


def normalize_dict(y_dict, chroms):
    robust_min = y_dict['robust_min']
    robust_max = y_dict['robust_max']

    y_dict_norm = {}
    for c in chroms:
        y_dict_norm[c] = (y_dict[c] - robust_min) / robust_max

    return y_dict_norm


def score(y_pred_dict, y_true_dict, chroms,
          gene_annotations, enh_annotations,
          window_size=25, prom_loc=80,
          y_var_true_dict=None):
    """Calculate score

    Args:
        bootstrapped_label_dict:
            { chr: label }
            Dict of labels of a bootstrapped (subsampled) sample
            For example of the full original label = [0,1,2,3] with
            4-fold bootstrapping without shuffling.
            Bootstrapped labels are like the following:
                [0,1,2], [0,1,3], [0,2,3], [1,2,3]
            We need to score a submission for each bootstrap index
            In real score calculation, we shuffle it with fixed
            random seed.
    """
    # concat all chromosomes
    y_pred_dict_norm = normalize_dict(y_pred_dict, chroms)
    y_true_dict_norm = normalize_dict(y_true_dict, chroms)

    y_pred = dict_to_arr(y_pred_dict, chroms)
    y_true = dict_to_arr(y_true_dict, chroms)

    y_pred_norm = dict_to_arr(y_pred_dict_norm, chroms)
    y_true_norm = dict_to_arr(y_true_dict_norm, chroms)

    if y_var_true_dict is None:
        y_var_true = None
    else:
        y_var_true = dict_to_arr(y_var_true_dict, chroms)

    output = Score(
        mse=mse(y_true_norm, y_pred_norm),
        gwcorr=gwcorr(y_true, y_pred),
        gwspear=gwspear(y_true, y_pred),
        mseprom=mseprom(y_true_dict_norm, y_pred_dict_norm, chroms,
                        gene_annotations,
                        window_size, prom_loc),
        msegene=msegene(y_true_dict_norm, y_pred_dict_norm, chroms,
                        gene_annotations,
                        window_size),
        mseenh=mseenh(y_true_dict_norm, y_pred_dict_norm, chroms,
                      enh_annotations,
                      window_size),
        msevar=msevar(y_true_norm, y_pred_norm, var=y_var_true),
        mse1obs=mse1obs(y_true_norm, y_pred_norm),
        mse1imp=mse1imp(y_true_norm, y_pred_norm),
    )
    return output


def write_to_db(score_db_record, db_file):
    cols = []
    values = []
    for attr in score_db_record._fields:
        cols.append(str(attr))
        val = getattr(score_db_record, attr)
        if isinstance(val, str):
            val = '"' + val + '"'
        else:
            val = str(val)
        values.append(val)

    query = DB_QUERY_INSERT.format(
        table=DB_TABLE_SCORE, cols=','.join(cols), values=','.join(values))
    log.info('SQL query: {}'.format(query))
    while True:
        try:
            conn = sqlite3.connect(db_file)
            c = conn.cursor()
            c.execute(query)
            c.close()
            conn.commit()
            conn.close()
        except sqlite3.OperationalError as e:
            print(e)
            conn.close()
            time.sleep(1)
            continue
        else:
            break    


def parse_arguments():
    import os
    py_path = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser(
        prog='ENCODE Imputation Challenge scoring script. .npy or .npz must be built '
             'with a correct --window-size. i.e. containing a value for each bin')
    parser.add_argument('npy_pred',
                        help='Submission .npy file to be scored')
    parser.add_argument('npy_true',
                        help='Truth .npy file')
    parser.add_argument('--var-npy',
                        help='Truth .npy file filled with a variance for each bin '
                        'instead of a raw signal value. '
                        'Variance must be computed over all cell types for a target '
                        'assay type. This should have a object form of '
                        '{ chrom: [var1, var2, ...] } length([var, ...]) should match '
                        'with the number of bins per chromosome. '
                        'Use build_var_npy.py TRUTH_NPY_CELL1 TRUTH_NPY_CELL2 ... '
                        'to build such .npy over all cell types for a target assay')
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
    p_score.add_argument('--bootstrap-chrom', nargs='*', default=[],
                         help='Bootstrapped chromosome groups. '
                              'Delimiter is whitespace for groups and '
                              'comma(,) in each group. Order is important.'
                              'e.g. "chr1,chr2 chr1,chrX chr2,chrX" means '
                              'three groups: (chr1,chr2), (chr1,chrX), (chr2,chrX)')
    p_score.add_argument('--gene-annotations',
                         default=os.path.join(
                            py_path,
                            'annot/hg38/gencode.v29.genes.gtf.bed.gz'),
                         help='Gene annotations BED file')
    p_score.add_argument('--enh-annotations',
                         default=os.path.join(
                            py_path,
                            'annot/hg38/F5.hg38.enhancers.bed.gz'),
                         help='Enhancer annotations BED file ')
    p_score.add_argument('--blacklist-file',
                         default=os.path.join(
                            py_path,
                            'annot/hg38/hg38.blacklist.bed.gz'),
                         help='Blacklist BED file. Bootstrap label will be '
                              'generated after removing overlapping regions '
                              'defined in this file.')
    p_score.add_argument('--window-size', default=25, type=int,
                         help='Window size for bigwig in bp')
    p_score.add_argument('--prom-loc', default=80, type=int,
                         help='Promoter location in a unit of window size '
                              '(--window-size). This is not in bp')
    p_score.add_argument('--validated', action='store_true',
                         help='For validated submissions (not for truth bigwigs)'
                              'with fixed interval length of 25 and valid '
                              'chromosome lengths. It will skip interpolation. '
                              'For truth bigwigs, it is recommended to convert '
                              'them into npy\'s or npz\'s by using '
                              'build_npy_from_bigwig.py')
    p_score.add_argument('--normalize-with-robust-min-max', action='store_true',
                         help='Normalize with robust min max.')
    p_out = parser.add_argument_group(
                        title='Output to file (TSV or DB)')
    p_out.add_argument('--out-file', default='output.tsv',
                       help='Write scores to TSV file')
    p_out.add_argument('--out-db-file',
                       help='Write metadata/scores to SQLite DB file')
    p_meta = parser.add_argument_group(
                        title='Submission metadata, which will be written to '
                              'DB together with scores. This will be used for '
                              'ranking later')
    p_meta.add_argument('--cell', '-c',
                        help='Cell type. e.g. C01')
    p_meta.add_argument('--assay', '-a',
                        help='Assay or histone mark. e.g. M01')
    p_meta.add_argument('--team-id', '-t', type=int,
                        help='Team ID (unique ID from Synapse)')
    p_meta.add_argument('--submission-id', '-s', type=int,
                        help='Submission ID (unique ID from Synapse)')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING',
                                 'CRITICAL', 'ERROR', 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    # some submission files have whitespace in path...
    if args.npy_pred is not None:
        args.npy_pred = args.npy_pred.strip("'")
    if args.npy_true is not None:
        args.npy_true = args.npy_true.strip("'")
    if args.var_npy is not None:
        args.var_npy = args.var_npy.strip("'")
    if args.out_file is not None:
        args.out_file = args.out_file.strip("'")

    if args.out_db_file is not None:
        args.out_db_file = args.out_db_file.strip("'")
        if not os.path.exists(args.out_db_file):
            raise ValueError('DB file does not exists')
        # if args.chrom != ['all']:
        #     raise ValueError(
        #         'Define "--chrom all" to write scores to DB file')
        if args.cell is None:
            raise ValueError(
                'Define "--cell CXX" to write scores to DB file')
        if args.assay is None:
            raise ValueError(
                'Define "--assay MXX" to write scores to DB file')
        if args.team_id is None:
            raise ValueError(
                'Define "--team-id" to write scores to DB file')
        if args.submission_id is None:
            raise ValueError(
                'Define "--submission-id" to write scores to DB file')

    if args.chrom == ['all']:
        args.chrom = ['chr' + str(i) for i in range(1, 23)] + ['chrX']
    args.chrom = sorted(args.chrom)

    if len(args.bootstrap_chrom) == 0:
        args.bootstrap_chrom = [(-1, args.chrom)]
    else:
        for i, _ in enumerate(args.bootstrap_chrom):
            args.bootstrap_chrom[i] = (i, args.bootstrap_chrom[i].split(','))
    print(args.bootstrap_chrom)

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def main():
    # read params
    args = parse_arguments()

    gc.disable()

    y_pred_dict = build_npy_from_bigwig(args.npy_pred, args.chrom,
                                        args.window_size, args.blacklist_file,
                                        args.validated)

    y_true_dict = build_npy_from_bigwig(args.npy_true, args.chrom,
                                        args.window_size, args.blacklist_file)

    if args.var_npy is None:
        y_var_true_dict = None
    elif args.var_npy.endswith(('.npy', '.npz')):
        log.info('Opening truth var numpy array file...')
        y_var_true_dict = numpy.load(args.var_npy, allow_pickle=True)[()]
    else:
        raise ValueError('Var true file should be a binned .npy or .npz.')

    log.info('Reading from enh_annotations...')
    enh_annotations = read_annotation_bed(args.enh_annotations)

    log.info('Reading from gene_annotations...')
    gene_annotations = read_annotation_bed(args.gene_annotations)

    with open(args.out_file, 'w') as fp:
        for k, bootstrap_chrom in args.bootstrap_chrom:
            log.info('Calculating score for bootstrap {} case...'.format(k))

            score_output = score(y_pred_dict, y_true_dict, bootstrap_chrom,
                           gene_annotations, enh_annotations,
                           args.window_size, args.prom_loc,
                           y_var_true_dict)
            # write to TSV
            s = "\t".join(['bootstrap_'+str(k)]+[str(o) for o in score_output])
            fp.write(s+'\n')
            print(s)

            # write to DB
            if args.out_db_file is not None:
                score_db_record = ScoreDBRecord(
                    args.submission_id,
                    args.team_id,
                    os.path.basename(args.npy_pred),
                    args.cell,
                    args.assay,
                    k,
                    *score_output)
                write_to_db(score_db_record, args.out_db_file)
            gc.collect()

    log.info('All done')


if __name__ == '__main__':
    main()

