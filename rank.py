#!/usr/bin/env python3
"""Imputation challenge rankigg script

Author:
    Jin Lee (leepc12@gmail.com)    
"""

import sys
import argparse
import sqlite3
import numpy
import time
from collections import namedtuple, defaultdict
from score import DB_TABLE_SCORE, ScoreDBRecord
from scipy.stats import rankdata
import logging

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)


DB_QUERY_GET = 'SELECT * FROM {table} ORDER BY bootstrap_id, submission_id;'
CELL_NAME = {
    'C02': 'adrenal_gland',
    'C20': 'heart_left_ventricle',
    'C35': 'omental_fat_pad',
    'C45': 'testis',
    'C05': 'BE2C',
    'C06': 'brain_microvascular_endothelial_cell',
    'C07': 'Caco-2',
    'C08': 'cardiac_fibroblast',
    'C11': 'dermis_microvascular_lymphatic_vessel_endothelial_cell',
    'C15': 'G401',
    'C16': 'GM06990',
    'C17': 'H1-hESC',
    'C18': 'H9',
    'C19': 'HAP-1',
    'C21': 'hematopoietic_multipotent_progenitor_cell',
    'C22': 'HL-60',
    'C23': 'IMR-90',
    'C24': 'K562',
    'C27': 'mesenchymal_stem_cell',
    'C28': 'MG63',
    'C30': 'NCI-H460',
    'C32': 'neural_stem_progenitor_cell',
    'C33': 'occipital_lobe',
    'C39': 'SJCRH30',
    'C40': 'SJSA1',
    'C41': 'SK-MEL-5',
    'C42': 'skin_fibroblast',
    'C43': 'skin_of_body',
    'C44': 'T47D',
    'C46': 'trophoblast_cell',
    'C47': 'upper_lobe_of_left_lung',
    'C48': 'urinary_bladder',
    'C49': 'uterus',
    'C51': 'WERI-Rb-1',
    'C12': 'DND-41',
    'C25': 'KMS-11',
    'C31': 'NCI-H929',
    'C34': 'OCI-LY7',
    'C01': 'adipose_tissue',
    'C03': 'adrenal_gland_embryonic',
    'C09': 'CD4-positive_alpha-beta_memory_T_cell',
    'C10': 'chorion',
    'C13': 'endocrine_pancreas',
    'C36': 'peripheral_blood_mononuclear_cell',
    'C37': 'prostate',
    'C38': 'RWPE2',
    'C50': 'vagina',
    'C04': 'amnion',
    'C29': 'myoepithelial_cell_of_mammary_gland',
    'C14': 'ES-I3',
    'C26': 'lower_leg_skin'
}
ASSAY_NAME = {
    'M01': 'ATAC-seq',
    'M02': 'DNase-seq',
    'M03': 'H2AFZ',
    'M04': 'H2AK5ac',
    'M05': 'H2AK9ac',
    'M06': 'H2BK120ac',
    'M07': 'H2BK12ac',
    'M08': 'H2BK15ac',
    'M09': 'H2BK20ac',
    'M10': 'H2BK5ac',
    'M11': 'H3F3A',
    'M12': 'H3K14ac',
    'M13': 'H3K18ac',
    'M14': 'H3K23ac',
    'M15': 'H3K23me2',
    'M16': 'H3K27ac',
    'M17': 'H3K27me3',
    'M18': 'H3K36me3',
    'M19': 'H3K4ac',
    'M20': 'H3K4me1',
    'M21': 'H3K4me2',
    'M22': 'H3K4me3',
    'M23': 'H3K56ac',
    'M24': 'H3K79me1',
    'M25': 'H3K79me2',
    'M26': 'H3K9ac',
    'M27': 'H3K9me1',
    'M28': 'H3K9me2',
    'M29': 'H3K9me3',
    'M30': 'H3T11ph',
    'M31': 'H4K12ac',
    'M32': 'H4K20me1',
    'M33': 'H4K5ac',
    'M34': 'H4K8ac',
    'M35': 'H4K91ac'
}

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

GlobalScore = namedtuple(
    'GlobalScore',
    ['team_id', 'name', 'score_lb', 'score_mean', 'score_ub', 'rank'])



TEAM_ID_FOR_ROUND1 = {
    0: 'Avocado_p0',
    1: 'Avocado_p1',
    2: 'Avocado_p2',
    3: 'Avocado_p3',
    4: 'Avocado_p4',
    5: 'Avocado_p5',
    6: 'Avocado_p6',
    7: 'Avocado_p7',
    8: 'Avocado_p8',
    9: 'Avocado_p9',
    10: 'Avocado_p10',
    20: 'Average',
    1000 : 'BrokenNodes',
    1010 : 'ENCODE_DCC_Imputes',
    1020 : 'Hongyang_Li_and_Yuanfang_Guan',
    1030 : 'KKT-ENCODE-Impute-model_1',
    1031 : 'KKT-ENCODE-Impute-model_2',
    1040 : 'LiPingChun',
    1050 : 'MLeipzig',
    1060 : 'noml'
}


def get_team_name(team_id):
    return TEAM_ID_FOR_ROUND1[team_id]


def get_cell_name(cell_id):
    if cell_id in CELL_NAME:
        return CELL_NAME[cell_id].replace('_', ' ')
    else:
        return str(cell_id)


def get_assay_name(assay_id):
    if assay_id in ASSAY_NAME:
        return ASSAY_NAME[assay_id].replace('_', ' ')
    else:
        return str(assay_id)


def score_record_factory(cursor, row):
    return ScoreDBRecord(*row)


def read_scores_from_db(db_file, chroms):
    """Read all rows by matching chromosomes
    Args:
        chroms:
            List of chromosome used for scoring. This will be
            converted into a comma-separated string and only
            rows with matching "chroms" field will be retrieved. 
            This is to filter out records scored with different
            set of chromosomes.
    Returns:
        All rows ordered by bootstrap_id and team_id
    """
    # valid_chrs_str = ','.join(sorted(chroms))
    query = DB_QUERY_GET.format(table=DB_TABLE_SCORE)  #, chroms=valid_chrs_str)
    log.info(query)

    while True:
        try:
            conn = sqlite3.connect(db_file)
            conn.row_factory = score_record_factory
            c = conn.cursor()
            c.execute(query)
            result = c.fetchall()
            c.close()
            conn.close()
        except sqlite3.OperationalError as e:
            print(e)
            conn.close()
            time.sleep(1)
            continue
        else:
            break
    
    if len(result) == 0:
        print('No records found. '
            'Did you forget to specify "--chrom"?')

    return result


def calc_combined_ranks(rows, measures_to_use):
    """Calculate ranks for combined measures
    """
    # make sure that the team_id is unique
    team_ids = [x.team_id for x in rows]
    submission_ids = [x.submission_id for x in rows]
    # submission_ids = [0 for x in rows]

    scores = numpy.zeros(len(team_ids), dtype=float)
    for user_i, attr in enumerate(measures_to_use):
        rank_method = RANK_METHOD_FOR_EACH_METRIC[attr]
        attr_scores = numpy.array([getattr(x, attr) for x in rows])
        if rank_method == 'ASCENDING':
            ranks = rankdata(-attr_scores, "average")
        elif rank_method == 'DESCENDING':
            ranks = rankdata(attr_scores, "average")
        else:
            raise Exception('Unknown rank_method.')

        pval_scores = numpy.log(ranks/float(len(ranks) + 1))
        scores += pval_scores
    ranks = rankdata(scores, "average")

    return dict(zip(zip(team_ids, submission_ids), ranks))


def calc_global_ranks(rows, measures_to_use):
    """Calculate global ranks

    Outputs:
        Markdown table for ranks
    """

    sample_grpd_results = defaultdict(lambda: defaultdict(list))
    all_users = set()
    for x in rows:
        sample_key = (x.cell, x.assay)
        sample_grpd_results[(x.cell, x.assay)][x.bootstrap_id].append(x)
        all_users.add(x.team_id)

    # group all submissions by cell and assay
    rv = {}
    global_scores = defaultdict(lambda: defaultdict(list))
    for (cell, assay), bootstrapped_submissions in sample_grpd_results.items():
        user_ranks = defaultdict(list)
        for index, submissions in bootstrapped_submissions.items():
            ranks = calc_combined_ranks(submissions, measures_to_use)

            obs_users = set(x[0] for x in ranks.keys())
            for (team_id, submission_id), rank in ranks.items():
                # print team_id, rank
                user_ranks[(team_id, submission_id)].append(rank)
                # user_ranks[(team_id, 0)].append(rank)
                global_scores[index][team_id].append(
                    min(0.5, rank/(len(ranks)+1))
                )
            for team_id in all_users - obs_users:
                global_scores[index][team_id].append(0.5)                

        print('# {} {} ({} {})'.format(cell, get_cell_name(cell), assay, get_assay_name(assay)))
        print(' | '.join(('Team', 'name', 'rank')))
        print('|'.join(('----',)*3))
        for (team_id, submission_id), ranks in sorted(
                user_ranks.items(), key=lambda x: sorted(x[1])[1]):
            print('%d | %s | %.2f' % (team_id, get_team_name(team_id), sorted(ranks)[1]))
        print()

    # group the scores by user
    user_grpd_global_scores = defaultdict(list)
    user_grpd_global_ranks = defaultdict(list)
    for bootstrap_id, bootstrap_global_scores in global_scores.items():
        sorted_scores = sorted(
            bootstrap_global_scores.items(), key=lambda x: sum(x[1]))
        ranks = rankdata([sum(x[1]) for x in sorted_scores])
        for (team_id, scores), rank in zip(sorted_scores, ranks):
            user_grpd_global_scores[team_id].append(sum(scores)/float(len(scores)))
            user_grpd_global_ranks[team_id].append(rank)
    global_data = []
    for team_id, scores in sorted(
            user_grpd_global_scores.items(), key=lambda x: sum(x[1])):
        global_data.append(GlobalScore(*[
            team_id, get_team_name(team_id), 
            min(scores), sum(scores)/len(scores), max(scores), 
            sorted(user_grpd_global_ranks[team_id])[1]
        ]))
    global_data = sorted(global_data, key=lambda x: (x.rank, x.score_mean))

    print()
    print('# Overall Results')
    print(' | '.join(('Team name', 'rank', 'Lower bound',
                      'Mean', 'Upperbound')))
    print('|'.join(('----',)*6))
    for x in global_data: 
        print('%s | %.2f | %.2f | %.2f | %.2f' % (
            x.name, x.rank, x.score_lb, x.score_mean, x.score_ub))

    return rv, global_data


def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE Imputation Challenge ranking'
                                          'script.')
    parser.add_argument('db_file',
                        help='SQLite3 DB file with all scores.')
    parser.add_argument('--show-score-only', action='store_true',
                        help='Show score (from DB) only')
    parser.add_argument('--chrom', nargs='+',
                        default=['all'],
                        help='List of chromosomes to be used for ranking')
    parser.add_argument('--measures-to-use', nargs='+',
                        default=['mse', 'mse1obs', 'mse1imp', 'gwcorr',
                                 'mseprom', 'msegene', 'mseenh'],
                        help='List of performance measures to be used for ranking')
    args = parser.parse_args()

    if args.chrom == ['all']:
        args.chrom = ['chr' + str(i) for i in range(1, 23)] + ['chrX']

    args.log_level = 'CRITICAL'

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def show_score(rows):
    print('\t'.join(['submission_id', 'team', 'cell_id', 'cell', 'assay_id', 'assay', 'bootstraip_id', 
                    'mse', 'gwcorr', 'gwspear', 'mseprom', 'msegene', 'mseenh',
                    'msevar', 'mse1obs', 'mse1imp']))
    for x in rows:
        mse = x.mse
        gwcorr = x.gwcorr
        gwspear = x.gwspear
        mseprom = x.mseprom
        msegene = x.msegene
        mseenh = x.mseenh
        msevar = x.msevar
        mse1obs = x.mse1obs
        mse1imp = x.mse1imp
        
        submission_id= x.submission_id
        team= get_team_name(x.team_id)
        cell_id = x.cell
        cell= get_cell_name(x.cell)
        assay_id = x.assay
        assay= get_assay_name(x.assay)
        bootstrap_id= x.bootstrap_id

        print('\t'.join([str(i) for i in \
                             [submission_id, team, cell_id, cell, assay_id, assay, bootstrap_id,
                              mse, gwcorr, gwspear, mseprom, msegene, mseenh,
                              msevar, mse1obs, mse1imp]]))


def main():
    # read params
    args = parse_arguments()

    log.info('Reading from DB file...')
    rows = read_scores_from_db(args.db_file, args.chrom)

    if args.show_score_only:
        log.info('List all scores...')
        show_score(rows)
    else:
        log.info('Calculate ranks...')
        rv, global_data = calc_global_ranks(rows, args.measures_to_use)

    log.info('All done.')


if __name__ == '__main__':
    main()
