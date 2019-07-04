#!/usr/bin/env python3
"""Imputation challenge rankigg script

Author:
    Jin Lee (leepc12@gmail.com)    
"""

import numpy
from collections import namedtuple, defaultdict
from scipy.stats import rankdata
from logger import log
from score_metric import RANK_METHOD_FOR_EACH_METRIC
from db import DB_QUERY_GET, read_scores_from_db
from challenge_metadata import get_cell_name, get_assay_name
from synapse import get_team_name


GlobalScore = namedtuple(
    'GlobalScore',
    ['team_id', 'name', 'score_lb', 'score_mean', 'score_ub', 'rank'])


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


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(
        description='ENCODE Imputation Challenge ranking script.')
    parser.add_argument('db_file',
                        help='SQLite3 DB file with all scores.')
    parser.add_argument('--show-score-only', action='store_true',
                        help='Show score (from DB) only')
    parser.add_argument('--chrom', nargs='+',
                        default=['all'],
                        help='List of chromosomes to be used for ranking')
    parser.add_argument('--measures-to-use', nargs='+',
                        default=['mse', 'gwcorr', 'gwspear', 'mseprom',
                                 'msegene', 'mseenh', 'msevar', 'mse1obs',
                                 'mse1imp'],
                        help='List of performance measures to be used for ranking')
    args = parser.parse_args()

    if args.chrom == ['all']:
        args.chrom = ['chr' + str(i) for i in range(1, 23)] + ['chrX']

    return args


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
