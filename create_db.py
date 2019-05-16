#!/usr/bin/env python3
"""Imputation challenge SQLite3 DB creation script
Author:
    Jin Lee (leepc12@gmail.com)
"""

import os
import sys
import argparse
import sqlite3
import logging
import score

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)


SCORE_DB_RECORD_VAR_TYPE = score.ScoreDBRecord(    
    submission_id='integer NOT NULL',
    team_id='integer NOT NULL',
    submission_fname='text NOT NULL',
    cell='text NOT NULL',
    assay='text NOT NULL',
    bootstrap_id='integer NOT NULL',
    chroms='text NOT NULL',
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


def parse_arguments():
    import os
    py_path = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser(
        prog='ENCODE Imputation Challenge SQLite3 Database creator.')
    parser.add_argument('db_file',
                        help='DB file.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING',
                                 'CRITICAL', 'ERROR', 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    if os.path.exists(args.db_file):
        raise ValueError('DB file already exists.')

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def main():
    # read params
    args = parse_arguments()

    log.info('Creating database...')

    try:
        conn = sqlite3.connect(args.db_file)
        c = conn.cursor()
        c.execute('CREATE TABLE IF NOT EXISTS {} ({});'.format(
            score.DB_TABLE_SCORE,
            ','.join([attr + ' ' + getattr(SCORE_DB_RECORD_VAR_TYPE, attr)
                        for attr in SCORE_DB_RECORD_VAR_TYPE._fields])))
    except Exception as e:
        print(e)
        sys.exit(1)
    finally:
        conn.close()

    log.info('All done.')


if __name__ == '__main__':
    main()
