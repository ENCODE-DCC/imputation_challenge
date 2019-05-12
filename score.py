# Imputed track evaluations
# Author: Jacob Schreiber, Jin Lee
# Contact: jmschreiber91@gmail.com, leepc12@gmail.com

import sys
import argparse
import numpy
import pyBigWig
import gzip
from collections import namedtuple
from sklearn.metrics import roc_auc_score
from scipy.stats import norm
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


def mseprom(y_true_dict, y_pred_dict, chroms,
            gene_annotations,
            window_size=25, prom_loc=80,
            bootstrapped_label_dict=None):
    """
    Args:
        y_true_dict: truth vector per chromosome
            { chr: y_true } where y_true is a numpy 1-dim array.

        y_pre_dict: predicted vector per chromosome
            { chr: y_pred } where y_pred is a numpy 1-dim array.

        bootstrapped_label_dict: bootstrapped_label per chromosome
    """
    sse, n = 0., 0.

    for chrom in chroms:
        y_true = y_true_dict[chrom]
        y_pred = y_pred_dict[chrom]
        if bootstrapped_label_dict is None:
            bootstrapped_label = None
        else:
            bootstrapped_label = set(bootstrapped_label_dict[chrom])

        for line in gene_annotations:
            chrom_, start, end, _, _, strand = line.split()
            start = int(start) // window_size
            end = int(end) // window_size + 1

            # if chrom_ in ('chrX', 'chrY', 'chrM'):
            #     continue

            if chrom_ != chrom:
                continue

            if bootstrapped_label is None:                
                if strand == '+':
                    sse += ((y_true[start-prom_loc: start] -
                             y_pred[start-prom_loc: start]) ** 2).sum()
                    n += y_true[start-prom_loc: start].shape[0]

                else:
                    sse += ((y_true[end: end+prom_loc] -
                             y_pred[end: end+prom_loc]) ** 2).sum()
                    n += y_true[end: end+prom_loc].shape[0]
            else:
                if strand == '+':
                    s = start-prom_loc
                    e = start
                else:
                    s = end
                    e = end+prom_loc
                for i in range(s, e):
                    if i not in bootstrapped_label:
                        continue
                    x = (y_true[i] - y_pred[i])
                    sse += (x*x)
                    n += 1

    return sse / n


def msegene(y_true_dict, y_pred_dict, chroms,
            gene_annotations,
            window_size=25, bootstrapped_label_dict=None):
    sse, n = 0., 0.

    for chrom in chroms:
        y_true = y_true_dict[chrom]
        y_pred = y_pred_dict[chrom]
        if bootstrapped_label_dict is None:
            bootstrapped_label = None
        else:
            bootstrapped_label = set(bootstrapped_label_dict[chrom])

        for line in gene_annotations:
            chrom_, start, end, _, _, strand = line.split()
            start = int(start) // window_size
            end = int(end) // window_size + 1

            # if chrom_ in ('chrX', 'chrY', 'chrM'):
            #     continue

            if chrom_ != chrom:
                continue
            if bootstrapped_label is None:
                sse += ((y_true[start:end] - y_pred[start:end]) ** 2).sum()
                n += end - start
            else:
                for i in range(start, end):
                    if i not in bootstrapped_label:
                        continue
                    x = (y_true[i] - y_pred[i])
                    sse += (x*x)
                    n += 1

    return sse / n


def mseenh(y_true_dict, y_pred_dict, chroms,
           enh_annotations,
           window_size=25, bootstrapped_label_dict=None):
    sse, n = 0., 0.

    for chrom in chroms:
        y_true = y_true_dict[chrom]
        y_pred = y_pred_dict[chrom]
        if bootstrapped_label_dict is None:
            bootstrapped_label = None
        else:
            bootstrapped_label = set(bootstrapped_label_dict[chrom])

        for line in enh_annotations:
            chrom_, start, end, _, _, _, _, _, _, _, _, _ = line.split()
            start = int(start) // window_size
            end = int(end) // window_size + 1

            # if chrom_ in ('chrX', 'chrY', 'chrM'):
            #     continue

            if chrom_ != chrom:
                continue
            if bootstrapped_label is None:
                sse += ((y_true[start:end] - y_pred[start:end]) ** 2).sum()
                n += end - start
            else:
                for i in range(start, end):
                    if i not in bootstrapped_label:
                        continue
                    x = (y_true[i] - y_pred[i])
                    sse += (x*x)
                    n += 1

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
        The mean-squared error that's weighted at each position by the
        variance.
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
    idealized vector of test cell type specific activation [1,0,0,0,....] as
    well as an idealized vector for test cell-type specific repression
    [0,1,1,1,.....].
    We use as the weight the maximum of these two mutual information values.
    As in MSEvar, the weights are normalized across all bins in the genome to
    sum to 1. The squared error between the predicted and true value at each
    genomic position is multiplied by this weight before being averaged across
    the genome.

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
        The mean-squared error that's weighted at each position by the
        variance.
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

    return (-norm.logpdf(y_pred, mu, std)).dot((y_true -
                                                y_pred) ** 2) / n


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


Score = namedtuple(
    'Score',
    ('mse', 'mse1obs', 'mse1imp', 'gwcorr', 'match1',
     'catch1obs', 'catch1imp', 'aucobs1', 'aucimp1',
     'mseprom', 'msegene', 'mseenh'))


ScoreDBRecord = namedtuple(
    'ScoreDBRecord',
    ('cell', 'assay', 'team_id', 'submission_id',
     'submission_fname', 'bootstrap_id') + Score._fields)


def score(y_pred_dict, y_true_dict, chroms,
          gene_annotations, enh_annotations,
          window_size=25, prom_loc=80,
          bootstrapped_label_dict=None):
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
    y_pred = dict_to_arr(y_pred_dict, bootstrapped_label_dict)
    y_true = dict_to_arr(y_true_dict, bootstrapped_label_dict)

    output = Score(
        mse=mse(y_true, y_pred),
        mse1obs=mse1obs(y_true, y_pred),
        mse1imp=mse1imp(y_true, y_pred),
        gwcorr=gwcorr(y_true, y_pred),
        match1=match1(y_true, y_pred),
        catch1obs=catch1obs(y_true, y_pred),
        catch1imp=catch1imp(y_true, y_pred),
        aucobs1=aucobs1(y_true, y_pred),
        aucimp1=aucimp1(y_true, y_pred),
        mseprom=mseprom(y_true_dict, y_pred_dict, chroms,
                        gene_annotations,
                        window_size, prom_loc,
                        bootstrapped_label_dict),
        msegene=msegene(y_true_dict, y_pred_dict, chroms,
                        gene_annotations,
                        window_size,
                        bootstrapped_label_dict),
        mseenh=mseenh(y_true_dict, y_pred_dict, chroms,
                      enh_annotations,
                      window_size,
                      bootstrapped_label_dict),
    )
    return output


def bw_to_dict(bw, chrs, window_size=25):
    """
    Returns:
        { chr: [] } where [] is a numpy 1-dim array
    """
    result = {}
    for c in chrs:
        result_per_chr = []
        chrom_len = bw.chroms()[c]
        for step in range((chrom_len-1)//window_size+1):
            start = step*window_size
            end = min((step+1)*window_size, chrom_len)
            stat = bw.stats(c, start, end, exact=True)
            # print(start,end,stat)
            if stat[0] is None:
                result_per_chr.append(0)
            else:
                result_per_chr.append(stat[0])
        result[c] = numpy.array(result_per_chr)
    return result


def dict_to_arr(d, bootstrapped_label_dict=None):
    """Concat vectors in d
    """
    result = []
    for c in d:
        if bootstrapped_label_dict is None:
            result.extend(d[c])
        else:
            label_for_chr = bootstrapped_label_dict[c]
            result.extend(d[c][label_for_chr])
    return numpy.array(result)


def parse_arguments():
    import os
    py_path = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser(
        prog='ENCODE Imputation Challenge scoring script.')
    parser.add_argument('bw_pred',
                        help='Submission bigwig file to be scored.')
    parser.add_argument('bw_true',
                        help='Truth bigwig file.')
    p_score = parser.add_argument_group(
                        title='Scoring parameters')
    p_score.add_argument('--chrom', nargs='+', required=True,
                         default=['all'],
                         help='List of chromosomes to be combined to be '
                              'scored. '
                              'Set as "all" (default) to score for all '
                              'chromosomes. '
                              '(e.g. "all" or "chr3 chr21") '
                              'It should be "all" to write scores to DB file.')
    p_score.add_argument('--bootstrapped-label-npy',
                         help='Bootstrapped label data .npy file. If given, '
                              'all folds '
                              'of bootstrapped samples will be scored')
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
    p_score.add_argument('--window-size', default=25, type=int,
                         help='Window size for bigwig in bp.')
    p_score.add_argument('--prom-loc', default=80, type=int,
                         help='Promoter location in a unit of window size '
                              '(--window-size). This is not in bp.')
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
    p_meta.add_argument('--mark', '-m',
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
    if args.bw_pred is not None:
        args.bw_pred = args.bw_pred.strip("'")
    if args.bw_true is not None:
        args.bw_true = args.bw_true.strip("'")
    if args.out_file is not None:
        args.out_file = args.out_file.strip("'")

    if args.out_db_file is not None:
        raise NotImplementedError('NO DB SUPPORT YET')

        args.out_db_file = args.out_db_file.strip("'")
        if args.chrom != ['all']:
            raise ValueError(
                'Define "--chrom all" to write scores to DB file.')
        if args.cell is None:
            raise ValueError(
                'Define "--cell CXX" to write scores to DB file.')
        if args.mark is None:
            raise ValueError(
                'Define "--mark MXX" to write scores to DB file.')
        if args.team_id is None:
            raise ValueError(
                'Define "--team-id" to write scores to DB file.')
        if args.submission_id is None:
            raise ValueError(
                'Define "--submission-id" to write scores to DB file.')

    if args.chrom == ['all']:
        args.chrom = ['chr' + str(i) for i in range(1, 23)] + ['chrX']

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def read_annotation_bed(bed):
    result = []
    if bed.endswith('gz'):
        with gzip.open(bed, 'r') as infile:
            for line in infile:
                result.append(line.decode("ascii"))
    else:
        with open(bed, 'r') as infile:
            for line in infile:
                result.append(line)
    return result


def main():
    # read params
    args = parse_arguments()

    log.info('Opening bigwig files...')
    bw_true = pyBigWig.open(args.bw_true.strip("'"))
    bw_pred = pyBigWig.open(args.bw_pred.strip("'"))

    log.info('Reading from enh_annotations...')
    enh_annotations = read_annotation_bed(args.enh_annotations)

    log.info('Reading from gene_annotations...')
    gene_annotations = read_annotation_bed(args.gene_annotations)

    if args.bootstrapped_label_npy is not None:
        log.info('Reading bootstrapped labels .npy...')
        bootstrapped_labels = numpy.load(args.bootstrapped_label_npy,
                                         allow_pickle=True)
    else:
        bootstrapped_labels = [None]

    with open(args.out_file, 'w') as fp:

        # dict: { chr: vector }
        log.info('Converting bigwig to numpy array...')
        y_pred_dict = bw_to_dict(bw_pred, args.chrom, args.window_size)
        y_true_dict = bw_to_dict(bw_true, args.chrom, args.window_size)

        for k, label in enumerate(bootstrapped_labels):

            if label is None:
                bootstrap = 'no_bootstrap'
            else:
                bootstrap = 'bootstrap-{}'.format(k+1)
            log.info('Calculating score for {} case...'.format(bootstrap))

            output = score(y_pred_dict, y_true_dict, args.chrom,
                           gene_annotations, enh_annotations,
                           args.window_size, args.prom_loc,
                           label)
            # write to TSV
            s = "\t".join([bootstrap]+[str(o) for o in output])
            fp.write(s+'\n')
            print(s)

            # write to DB


    log.info('All done.')


if __name__ == '__main__':
    main()

"""DEV NOTE:

def write_to_db():
    for x in range(0, timeout):
        try:
            with connection:
            connection.execute(sql)
        except:
            time.sleep(1)
            pass
        finally:
            break
    else:
        with connection:
            connection.execute(sql)
"""
