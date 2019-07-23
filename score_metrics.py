#!/usr/bin/env python3
"""Imputation challenge metrics

Author:
    Jacob Schreiber (jmschreiber91@gmail.com)
    Jin Lee (leepc12@gmail.com)
"""

import numpy
from collections import namedtuple
from sklearn.metrics import roc_auc_score
from scipy.stats import norm, spearmanr
from logger import log


Score = namedtuple(
    'Score',
    ('mse', 'gwcorr', 'gwspear', 'mseprom', 'msegene', 'mseenh',
     'msevar', 'mse1obs', 'mse1imp')
)

# Ascending (the bigger the better) or descending order for each metric
RANK_METHOD_FOR_EACH_METRIC = {
    'mse': 'DESCENDING',
    'gwcorr': 'ASCENDING',
    'gwspear': 'ASCENDING',
    'mseprom': 'DESCENDING',
    'msegene': 'DESCENDING',
    'mseenh': 'DESCENDING',
    'msevar': 'DESCENDING',
    'mse1obs': 'DESCENDING',
    'mse1imp': 'DESCENDING'
}


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


def find_robust_min_max(x, pct_thresh=0.05, top_bottom_bin_range=2000000):
    y = x[x > 0]
    idxs = numpy.argsort(y)
    abs_max = y[idxs[-1]]
    abs_min = y[idxs[0]]
    robust_max = y[idxs[-int(pct_thresh * top_bottom_bin_range)]]
    robust_min = y[idxs[int(pct_thresh * top_bottom_bin_range)]]
    log.info('Array length original, non-zero: {}, {}'.format(len(x), len(y)))
    log.info('Absolute min, max: {}, {}'.format(abs_min, abs_max))
    log.info('Robust min, max: {}, {}'.format(robust_min, robust_max))
    return robust_min, robust_max


def normalize_dict(y_dict, chroms):
    #robust_min = y_dict['robust_min']
    #robust_max = y_dict['robust_max']
    #y_dict_norm = {}
    #for c in chroms:
    #    y = numpy.array(y_dict[c])
    #    y[y <= robust_min] = robust_min
    #    y[y >= robust_max] = robust_max
    #    y_dict_norm[c] = (y - robust_min) / robust_max
    #return y_dict_norm
    return y_dict
