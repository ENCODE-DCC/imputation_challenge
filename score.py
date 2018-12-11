# Imputed track evaluations
# Author: Jacob Schreiber, Jin Lee
# Contact: jmschreiber91@gmail.com, leepc12@gmail.com

import sys
import argparse
import numpy
import pyBigWig
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

    with open(gene_annotations, 'r') as infile:
        for line in infile:
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

    with open(gene_annotations, 'r') as infile:
        for line in infile:
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

    with open(enh_annotations, 'r') as infile:
        for line in infile:
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

    with open(args.out,'w') as fp:
        for chrom in args.chrom:
            log.info('Scoring for chrom {}...'.format(chrom))
            y_predicted_ = bw_predicted.values(chrom, 0, bw_predicted.chroms()[chrom])
            y_true_ = bw_true.values(chrom, 0, bw_true.chroms()[chrom])
            y_predicted = numpy.array(y_predicted_)
            y_true = numpy.array(y_true_)

            output = [chrom]
            for func in mse, mse1obs, mse1imp, gwcorr, match1, catch1obs, catch1imp, aucobs1, aucimp1:
                output.append(func(y_true, y_predicted))
            output.append(mseprom(y_true, y_predicted, args.gene_annotations, args.window_size, args.prom_loc))
            output.append(msegene(y_true, y_predicted, args.gene_annotations, args.window_size))
            output.append(mseenh(y_true, y_predicted, args.enh_annotations, args.window_size))
            fp.write("\t".join(output))
            print("\t".join(output))

    log.info('All done.')

if __name__=='__main__':
    main()
