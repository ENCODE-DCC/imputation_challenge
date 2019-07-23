#!/usr/bin/env python3
"""Imputation challenge sqlite3 DB functions

Author:
    Jin Lee (leepc12@gmail.com)
"""

import sys
import sqlite3
from collections import namedtuple
from score_metrics import Score
from logger import log


ScoreDBRecord = namedtuple(
    'ScoreDBRecord',
    ('submission_id', 'team_id', 'submission_fname', 'cell', 'assay',
     'bootstrap_id') + Score._fields
)

DB_TABLE_SCORE = 'score'
DB_QUERY_INSERT = 'INSERT INTO {table} ({cols}) VALUES ({values});'
#DB_QUERY_GET = 'SELECT * FROM {table} ORDER BY bootstrap_id, submission_id;'
# to select latest submission
DB_QUERY_GET = 'SELECT t.* FROM {table} t \
INNER JOIN ( \
	SELECT team_id, cell, assay, bootstrap_id, max(submission_id) as MaxSID \
	FROM {table} \
	GROUP BY team_id, cell, assay, bootstrap_id \
) tm WHERE t.team_id = tm.team_id AND t.cell = tm.cell \
AND t.assay = tm.assay AND t.bootstrap_id = tm.bootstrap_id \
AND t.submission_id = tm.MaxSID \
ORDER BY t.bootstrap_id, t.submission_id;'

SCORE_DB_RECORD_VAR_TYPE = ScoreDBRecord(
    submission_id='integer NOT NULL',
    team_id='integer NOT NULL',
    submission_fname='text NOT NULL',
    cell='text NOT NULL',
    assay='text NOT NULL',
    bootstrap_id='integer NOT NULL',
    mse='double NOT NULL',
    gwcorr='double NOT NULL',
    gwspear='double NOT NULL',
    mseprom='double NOT NULL',
    msegene='double NOT NULL',
    mseenh='double NOT NULL',
    msevar='double NOT NULL',
    mse1obs='double NOT NULL',
    mse1imp='double NOT NULL'
)


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
    def score_record_factory(cursor, row):
        return ScoreDBRecord(*row)
    
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
    
    #if len(result) == 0:
    #    print('No records found. '
    #        'Did you forget to specify "--chrom"?')

    return result


def create_db(db_file):
    log.info('Creating database...')

    try:
        conn = sqlite3.connect(db_file)
        c = conn.cursor()
        c.execute('CREATE TABLE IF NOT EXISTS {} ({});'.format(
            DB_TABLE_SCORE,
            ','.join([attr + ' ' + getattr(SCORE_DB_RECORD_VAR_TYPE, attr)
                        for attr in SCORE_DB_RECORD_VAR_TYPE._fields])))
    except Exception as e:
        print(e)
        sys.exit(1)
    finally:
        conn.close()

    log.info('All done.')


def parse_arguments():
    import argparse
    import os

    parser = argparse.ArgumentParser(
        description='ENCODE Imputation Challenge SQLite3 Database creator.')
    parser.add_argument('db_file', help='DB file.')
    args = parser.parse_args()

    if os.path.exists(args.db_file):
        raise ValueError('DB file already exists.')

    return args


def main():
    args = parse_arguments()
    create_db(args.db_file)


if __name__ == '__main__':
    main()
