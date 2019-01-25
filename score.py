# Imputed track evaluations
# Author: Jacob Schreiber, Jin Lee
# Contact: jmschreiber91@gmail.com, leepc12@gmail.com

import sys
import argparse
import numpy
import math
import pyBigWig
import gzip
from sklearn.metrics import roc_auc_score
import logging

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)

BIG_INT = 99999999

def mse(y_true, y_pred):
    return ((y_true - y_pred) ** 2.).mean()

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

def gwcorr(y_true, y_pred):
    return numpy.corrcoef(y_true, y_pred)[0, 1]

def match1(y_true, y_pred):
    n = int(y_true.shape[0] * 0.01)
    y_true_sorted = numpy.sort(y_true)
    y_pred_sorted = numpy.sort(y_pred)
    
    y_true_top1 = y_true_sorted[-n]
    y_pred_top1 = y_pred_sorted[-n]

    y_true_top = y_true >= y_true_top1
    y_pred_top = y_pred >= y_pred_top1

    return (y_true_top & y_pred_top).sum()

def catch1obs(y_true, y_pred):
    n = int(y_true.shape[0] * 0.01)
    y_true_sorted = numpy.sort(y_true)
    y_pred_sorted = numpy.sort(y_pred)
    
    y_true_top1 = y_true_sorted[-n]
    y_pred_top1 = y_pred_sorted[-n*5]

    y_true_top = y_true >= y_true_top1
    y_pred_top = y_pred >= y_pred_top1

    return (y_true_top & y_pred_top).sum()

def catch1imp(y_true, y_pred):
    n = int(y_true.shape[0] * 0.01)
    y_true_sorted = numpy.sort(y_true)
    y_pred_sorted = numpy.sort(y_pred)
    
    y_true_top1 = y_true_sorted[-n*5]
    y_pred_top1 = y_pred_sorted[-n]

    y_true_top = y_true >= y_true_top1
    y_pred_top = y_pred >= y_pred_top1

    return (y_true_top & y_pred_top).sum()

def aucobs1(y_true, y_pred):
    n = int(y_true.shape[0] * 0.01)
    y_true_sorted = numpy.sort(y_true)

    y_true_top1 = y_true_sorted[-n]
    y_true_top = y_true >= y_true_top1

    return roc_auc_score(y_true_top, y_pred)

def aucimp1(y_true, y_pred):
    n = int(y_true.shape[0] * 0.01)
    y_pred_sorted = numpy.sort(y_pred)

    y_pred_top1 = y_pred_sorted[-n]
    y_pred_top = y_pred >= y_pred_top1

    return roc_auc_score(y_pred_top, y_true)

def mseprom(y_true, y_pred, chrom, gene_annotations, window_size=25, prom_loc=80):
    sse, n = 0., 0.

    for line in gene_annotations:
        chrom_, start, end, _, _, strand = line.split()
        start = int(start) // window_size
        end = int(end) // window_size + 1
        
        # if chrom_ in ('chrX', 'chrY', 'chrM'):
        #     continue

        if chrom_ != chrom:
            continue

        if strand == '+':
            sse += ((y_true[start-prom_loc: start] - y_pred[start-prom_loc: start]) ** 2).sum()
            n += y_true[start-prom_loc: start].shape[0]

        else:
            sse += ((y_true[end: end+prom_loc] - y_pred[end: end+prom_loc]) ** 2).sum()
            n += y_true[end: end+prom_loc].shape[0]

    return sse / n

def msegene(y_true, y_pred, chrom, gene_annotations, window_size=25):
    sse, n = 0., 0.

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

def mseenh(y_true, y_pred, chrom, enh_annotations, window_size=25):
    sse, n = 0., 0.

    for line in enh_annotations:
        chrom_, start, end, _, _, _, _, _, _, _, _, _ = line.split()
        start = int(start) // window_size
        end = int(end) // window_size + 1
        
        # if chrom_ in ('chrX', 'chrY', 'chrM'):
        #     continue

        if chrom_ != chrom:
            continue

        sse += ((y_true[start:end] - y_pred[start:end]) ** 2).sum()
        n += end - start

    return sse / n

def mseVar(y_true, y_pred, y_all):
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

    y_all: numpy.ndarray, shape=(n_positions, n_celltypes)
        The true signal from all the cell types to calculate the variance over.

    Returns
    -------
    mse: float
        The mean-squared error that's weighted at each position by the variance.
    """

    var = numpy.std(y_all, axis=0) ** 2
    var /= var.sum()
    return ((y_true - y_pred) ** 2).dot(var)

def mseSpec(y_true, y_pred, y_all):
    """Calculates the MSE weighted by the specificity of the signal.

    Anshul has not sent me this yet so I'm not sure what to put here.

    According to the wiki: Test values that are cell type-specific are 
    typically the hardest to impute; hence, we upweight such bins. A 
    weight is computed for each bin in the genome based on how specific 
    the signal is in the test cell type relative to training cell types. 
    Specifically, for each bin in the genome, we compute the mutual information 
    between the signal vector for the bin across test and training cell types 
    [S_t, S_r1, S_r2, S_r3, ....] (where S_t is the signal in the test cell 
    type and S_ri is the signal value in the i-th training cell type) and an 
    idealized vector of test cell type specific activation [1,0,0,0,....] as well 
    as an idealized vector for test cell-type specific repression [0,1,1,1,.....]. 
    We use as the weight the maximum of these two mutual information values. As 
    in MSEvar, the weights are normalized across all bins in the genome to sum to 1. 
    The squared error between the predicted and true value at each genomic position 
    is multiplied by this weight before being averaged across the genome.

    Parameters
    ----------
    y_true: numpy.ndarray, shape=(n_positions,)
        The true signal

    y_pred: numpy.ndarray, shape=(n_positions,)
        The predicted signal

    y_all: numpy.ndarray, shape=(n_positions, n_celltypes)
        The true signal from all the cell types to calculate the variance over.

    Returns
    -------
    mse: float
        The mean-squared error that's weighted at each position by the variance.
    """
    
    return 0

def mseSpec2(y_true, y_pred, y_all):
    """Calculates the MSE weighted by the specificity of the signal.

    This is my implementation of a metric that weights each position by how
    different it is from the training cell types. Essentially, a normal
    distribution is fit to the values in the training cell types at each
    position, and the -logp is calculated for the test set value. This will
    give high weights to values that are very different, and low weights
    to values that are similar.

    Parameters
    ----------
    y_true: numpy.ndarray, shape=(n_positions,)
        The true signal

    y_pred: numpy.ndarray, shape=(n_positions,)
        The predicted signal

    y_all: numpy.ndarray, shape=(n_positions, n_celltypes)
        The true signal from all the cell types to calculate the normal
        disribution over.

    Returns
    -------
    mse: float
        The mean-squared error that's weighted at each position by -logp.
    """

    mu = y_all.mean(axis=0)
    std = y_all.std(axis=0)
    std[std < 0.01] = 0.01
    n = y_true.shape[0]

    return (-norm.logpdf(y_pred, mu, std)).dot((y_true - y_pred) ** 2) / n

def mseDiff(y_true, y_pred):
    """Calculates the MSE over the difference in the tracks.

    This metric captures the smoothness of the signal. The difference between
    adjacent positions is first calculated for both the real signal and the
    imputed signal, and then the MSE is calculated over that.

    Parameters
    ----------
    y_true: numpy.ndarray, shape=(n_positions,)
        The true signal

    y_pred: numpy.ndarray, shape=(n_positions,)
        The predicted signal

    Returns
    -------
    mse: float
        The mean-squared error over the difference.
    """

    return ((numpy.diff(y_true) - numpy.diff(y_pred)) ** 2).mean()

def bw_to_arr(bw, chrom, window_size):
    result = []    
    chrom_len = bw.chroms()[chrom]
    for step in range((chrom_len-1)//window_size+1):
        start = step*window_size
        end = min((step+1)*window_size,chrom_len)
        stat = bw.stats(chrom, start, end, exact=True)
        # print(start,end,stat)
        if stat[0]==None:
            result.append(0)
        else:
            result.append(stat[0])
    return numpy.array(result)

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE Imputation Challenge scoring script.',
                                        description='')
    parser.add_argument('bw_predicted', type=str,
                        help='Bigwig file to be scored.')
    parser.add_argument('bw_true', type=str,
                        help='Bigwig file to be scored.')
    parser.add_argument('--gene-annotations', type=str, required=True,
                        help='Gene annotations bed file (not in a compressed gz format).')
    parser.add_argument('--enh-annotations', type=str, required=True,
                        help='Enhancer annotations bed file (not in a compressed gz format).')
    parser.add_argument('--chrom', nargs='+', type=str, required=True,
                        help='List of chromosomes to score for (e.g. chr3 chr21).')
    parser.add_argument('--window-size', default=25, type=int,
                        help='Window size for bigwig in bp.')
    parser.add_argument('--prom-loc', default=80, type=int,
                        help='Promoter location in a unit of window size (--window-size). This is not in bp.')
    parser.add_argument('--out', default='output.tsv', type=str,
                            help='Output TSV file. (chrom, mse, mse1obs, mse1imp, gwcorr, match1, catch1obs, '
                                'catch1imp, aucobs1, aucimp1, mseprom, msegene, mseenh)')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--log-level', default='INFO', 
                        choices=['NOTSET','DEBUG','INFO','WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

def main():
    # read params
    args = parse_arguments()

    log.info('Opening bigwig files...')
    bw_true = pyBigWig.open(args.bw_true)
    bw_predicted = pyBigWig.open(args.bw_predicted)

    log.info('Reading from enh_annotations...')
    enh_annotations=[]
    if args.enh_annotations.endswith('gz'):
        with gzip.open(args.enh_annotations, 'r') as infile:
            for line in infile:
                enh_annotations.append(line.decode("ascii"))        
    else:
        with open(args.enh_annotations, 'r') as infile:
            for line in infile:
                enh_annotations.append(line)

    log.info('Reading from gene_annotations...')
    gene_annotations=[]

    if args.gene_annotations.endswith('gz'):
        with gzip.open(args.gene_annotations, 'r') as infile:
            for line in infile:
                gene_annotations.append(line.decode("ascii"))
    else:
        with open(args.gene_annotations, 'r') as infile:
            for line in infile:
                gene_annotations.append(line)

    # print(len(enh_annotations))
    # print(len(gene_annotations))
    # print(enh_annotations[0:2])
    # print(gene_annotations[0:2])
    # print(enh_annotations[0].split())
    # print(gene_annotations[0].split())
    # sys.exit(1)

    with open(args.out,'w') as fp:
        for chrom in args.chrom:
            log.info('Scoring for chrom {}...'.format(chrom))
            y_true = bw_to_arr(bw_true, chrom, args.window_size)
            log.info('y_true_len: {}'.format(len(y_true)))
            y_predicted = bw_to_arr(bw_predicted, chrom, args.window_size)
            log.info('y_predicted_len: {}'.format(len(y_predicted)))

            output = [chrom]
            for func in mse, mse1obs, mse1imp, gwcorr, match1, catch1obs, catch1imp, aucobs1, aucimp1:
                output.append(func(y_true, y_predicted))
            output.append(mseprom(y_true, y_predicted, chrom, gene_annotations, args.window_size, args.prom_loc))
            output.append(msegene(y_true, y_predicted, chrom, gene_annotations, args.window_size))
            output.append(mseenh(y_true, y_predicted, chrom, enh_annotations, args.window_size))
            fp.write("\t".join([str(o) for o in output]))
            fp.write("\n")
            print("\t".join([str(o) for o in output]))

    log.info('All done.')

if __name__=='__main__':
    main()
