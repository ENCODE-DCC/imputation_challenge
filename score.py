#!/usr/bin/env python3
"""Imputation challenge scoring script
Author:
    Jin Lee (leepc12@gmail.com)
"""

import os
import time
import gc
from score_metrics import Score, normalize_dict
from score_metrics import mse, mseprom, msegene, mseenh, msevar, mse1obs, mse1imp
from score_metrics import gwcorr, gwspear
from db import write_to_db, ScoreDBRecord
from bw_to_npy import load_bed, load_npy, bw_to_dict, dict_to_arr
import logging
log = logging.getLogger(__name__)


def parse_submission_filename(bw_file):
    """Filename should be CXXMYY.bigwig or CXXMYY.bw
    """
    basename = os.path.basename(bw_file)
    cell = basename[0:3]
    assay = basename[3:6]
    return cell, assay


def score(y_pred_dict, y_true_dict, chroms,
          gene_annotations, enh_annotations,
          window_size=25, prom_loc=80,
          y_var_dict=None):
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

    if y_var_dict is None:
        y_var_true = None
    else:
        y_var_true = dict_to_arr(y_var_dict, chroms)

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


def parse_arguments():
    import argparse
    import os

    py_path = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser(
        description='ENCODE Imputation Challenge scoring script. .npy or .npz must be built '
             'with a correct --window-size. i.e. containing a value for each bin')
    parser.add_argument('pred_npy_or_bw',
                        help='Submission .npy or .bigwig file to be scored')
    parser.add_argument('true_npy_or_bw',
                        help='Truth .npy or .bigwig file')
    parser.add_argument('--download-submissions-from-syn-eval-queue', action='store_true',
                         help='Download RECEIVED submissions from Synapse '
                              'evaluation queue.')
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
                         help='For validated submissions (not for truth bigwigs) '
                              'with fixed interval length of 25 and valid '
                              'chromosome lengths. It will skip interpolation. '
                              'For truth bigwigs, it is recommended to convert '
                              'them into npy\'s or npz\'s by using '
                              'bw_to_npy.py')
    #p_score.add_argument('--normalize-with-robust-min-max', action='store_true',
    #                     help='Normalize with robust min max.')
    p_out = parser.add_argument_group(
                        title='Output to file (TSV or DB)')
    p_out.add_argument('--db-file',
                       help='Write metadata/scores to SQLite DB file')
    p_meta = parser.add_argument_group(
                        title='Submission metadata, which will be written to '
                              'DB together with scores. This will be used for '
                              'ranking later')
    p_meta.add_argument('--team-id', '-t', type=int,
                        help='Team ID (unique ID from Synapse)')
    p_meta.add_argument('--submission-id', '-s', type=int,
                        help='Submission ID (unique ID from Synapse)')
    args = parser.parse_args()

    if args.db_file is not None:
        if not os.path.exists(args.db_file):
            raise ValueError('DB file does not exists')

    if args.chrom == ['all']:
        args.chrom = ['chr' + str(i) for i in range(1, 23)] + ['chrX']
    args.chrom = sorted(args.chrom)

    if len(args.bootstrap_chrom) == 0:
        args.bootstrap_chrom = [(-1, args.chrom)]
    else:
        for i, _ in enumerate(args.bootstrap_chrom):
            args.bootstrap_chrom[i] = (i, args.bootstrap_chrom[i].split(','))
    print(args.bootstrap_chrom)

    return args


def main():
    args = parse_arguments()

    y_pred_dict = bw_to_dict(args.pred_npy_or_bw, args.chrom,
                             args.window_size, args.blacklist_file)
    y_true_dict = bw_to_dict(args.true_npy_or_bw, args.chrom,
                             args.window_size, args.blacklist_file)
    if args.var_npy is None:
        y_var_dict = None
    elif args.var_npy.endswith(('.npy', '.npz')):
        log.info('Opening truth var numpy array file...')
        y_var_dict = load_npy(args.var_npy)
    else:
        raise ValueError('Var true file should be a binned .npy or .npz.')

    enh_annotations = load_bed(args.enh_annotations)
    gene_annotations = load_bed(args.gene_annotations)

    gc.disable()

    cell, assay = parse_submission_filename(args.pred_npy_or_bw)

    for k, bootstrap_chrom in args.bootstrap_chrom:
        log.info('Calculating score for bootstrap {} case...'.format(k))

        score_output = score(y_pred_dict, y_true_dict, bootstrap_chrom,
                             gene_annotations, enh_annotations,
                             args.window_size, args.prom_loc,
                             y_var_dict)
        s = "\t".join(['bootstrap_'+str(k)]+[str(o) for o in score_output])
        print(s)

        # write to DB
        if args.db_file is not None:
            score_db_record = ScoreDBRecord(
                args.submission_id,
                args.team_id,
                os.path.basename(bw),
                cell,
                assay,
                k,
                *score_output)
            write_to_db(score_db_record, args.db_file)
        gc.collect()

    log.info('All done')


if __name__ == '__main__':
    main()

