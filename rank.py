#!/usr/bin/env python3
"""Imputation challenge rankimg script

Author:
    Jin Lee (leepc12@gmail.com)
"""

import sys
import argparse
import sqlite3
import numpy
from collections import namedtuple, defaultdict
from score import DB_TABLE_SCORE, ScoreDBRecord
from scipy.stats import rankdata
import logging

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)


DB_QUERY_GET = 'SELECT * FROM {table} WHERE chroms="{chroms}" ORDER BY bootstrap_id, submission_id;'

ATTRS_TO_RANK = ('mse', 'mse1obs', 'mse1imp', 'gwcorr', 'match1',
                 'catch1obs', 'catch1imp', 'aucobs1', 'aucimp1',
                 'mseprom', 'msegene', 'mseenh')

GlobalScore = namedtuple(
    'GlobalScore',
    ['team_id', 'name', 'score_lb', 'score_mean', 'score_ub', 'rank'])



def get_team_name(team_id):
    return str(team_id)


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
    valid_chrs_str = ','.join(sorted(chroms))
    query = DB_QUERY_GET.format(table=DB_TABLE_SCORE, chroms=valid_chrs_str)
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
    return result


def calc_combined_ranks(rows):
    # make sure that the team_id is unique
    team_ids = [x.team_id for x in rows]
    submission_ids = [x.submission_id for x in rows]
    # submission_ids = [0 for x in rows]

    scores = numpy.zeros(len(team_ids), dtype=float)
    for user_i, attr in enumerate(ATTRS_TO_RANK):
        attr_scores = numpy.array([getattr(x, attr) for x in rows])
        ranks = rankdata(-attr_scores, "average")
        pval_scores = numpy.log(ranks/float(len(ranks) + 1))
        scores += pval_scores
    ranks = rankdata(scores, "average")

    return dict(zip(zip(team_ids, submission_ids), ranks))


def calc_global_ranks(rows):
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
            ranks = calc_combined_ranks(submissions)

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

    print('\t'.join(('name', 'rank', 'lb', 'mean', 'ub')))
    for x in global_data: 
        print('%s\t%.2f\t%.2f\t%.2f\t%.2f' % (
              x.name, x.rank, x.score_lb, x.score_mean, x.score_ub))
    print('# Overall Results')
    print(' | '.join(
         ('Team name', 'rank', 'Lower bound', 'Mean', 'Upperbound')))
    print('|'.join(('---',)*6))
    for x in global_data: 
        print('%s | %.2f | %.2f | %.2f | %.2f' % (
            x.name, x.rank, x.score_lb, x.score_mean, x.score_ub))

    return rv, global_data


def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE Imputation Challenge ranking'
                                          'script.')
    parser.add_argument('db_file',
                        help='SQLite3 DB file with all scores.')
    parser.add_argument('--chrom', nargs='+',
                        default=['all'],
                        help='List of chromosomes to be used for ranking')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING',
                                 'CRITICAL', 'ERROR', 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    if args.chrom == ['all']:
        args.chrom = ['chr' + str(i) for i in range(1, 23)] + ['chrX']

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def main():
    # read params
    args = parse_arguments()

    log.info('Reading from DB file...')
    rows = read_scores_from_db(args.db_file, args.chrom)

    log.info('Calculate ranks...')
    rv, global_data = calc_global_ranks(rows)

    log.info('All done.')


if __name__ == '__main__':
    main()
