#!/usr/bin/env python3
"""Imputation challenge validation script for round2
Author:
    Jin Lee (leepc12@gmail.com)
"""

import os
import re
import time
import shutil
import math
#import gc
import traceback
import synapseclient
import multiprocessing
from validate import validate
from io import StringIO
from score_leaderboard import mkdir_p, send_message
from score_leaderboard import WIKI_TEMPLATE_SUBMISSION_STATUS, RE_PATTERN_SUBMISSION_FNAME
from logger import log


BIG_INT = 99999999  # for multiprocessing


ROUND2_VALID_CELL_ASSAY = [
'C05M17',
'C05M18',
'C05M20',
'C05M29',
'C06M16',
'C06M17',
'C06M18',
'C07M20',
'C07M29',
'C12M01',
'C12M02',
'C19M16',
'C19M17',
'C19M18',
'C19M20',
'C19M22',
'C19M29',
'C22M16',
'C22M17',
'C28M17',
'C28M18',
'C28M22',
'C28M29',
'C31M01',
'C31M02',
'C31M16',
'C31M29',
'C38M01',
'C38M02',
'C38M17',
'C38M18',
'C38M20',
'C38M22',
'C38M29',
'C39M16',
'C39M17',
'C39M18',
'C39M20',
'C39M22',
'C39M29',
'C40M16',
'C40M17',
'C40M18',
'C40M20',
'C40M22',
'C40M29',
'C51M16',
'C51M17',
'C51M18',
'C51M20',
'C51M29'
]

def is_valid_round2_cell_assay(cell, assay):
    return (cell + assay) in ROUND2_VALID_CELL_ASSAY

def update_wiki_for_round2(syn, team_name_dict, args):
    # calculate ranks and update round2 wiki
    log.info('Updating wiki...')
    wiki_id_map = {
        k.split(':')[0]: k.split(':')[1] for k in args.round2_wiki_id.split(',')
    }
    for k, wiki_id in wiki_id_map.items():
        w = syn.getWiki(args.project_id, wiki_id)

        if k == 'submission_status':
            title = 'Submission status'
            markdown = WIKI_TEMPLATE_SUBMISSION_STATUS.format(
                eval_queue_id=args.eval_queue_id)

        else:
            raise Exception('invalid wiki type')

        w.markdown = markdown
        w.title = title
        w = syn.store(w)

    return None

def parse_arguments():
    import argparse
    import os

    py_path = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser(
        description='ENCODE Imputation Challenge round2 script. ')
    parser.add_argument('eval_queue_id',
                        help='Synapse evaluation queue ID to retreive submissions from.')
    p_score = parser.add_argument_group(
                        title='Scoring parameters')
    p_score.add_argument('--window-size', default=25, type=int,
                         help='Window size for bigwig in bp')
    p_score.add_argument('--update-wiki-only', action='store_true',
                         help='Update wiki based on DB file (--db-file) without '
                              'scoring submissions')
    p_sys = parser.add_argument_group(
                        title='System and resource settings')
    p_sys.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize scoring (per) submission')
    p_sys.add_argument('--team-name-tsv',
                        help='TSV file with team_id/team_name (1st col/2nd col).')
    p_syn = parser.add_argument_group(
                        title='Communitation with synapse')
    p_syn.add_argument('--dry-run', action='store_true',
                       help='Do not update submission\'s status on synapse.')
    p_syn.add_argument('--project-id', default='syn17083203',
                       help='Synapse project ID.')
    p_syn.add_argument('--round2-wiki-id',
                       default='submission_status:594309',
                       help='Comma-delimited Synapse wiki ID for round2.')
    p_syn.add_argument('--submission-dir', default='./submissions',
                       help='Download submissions here.')
    p_syn.add_argument('--send-msg-to-admin', action='store_true',
                       help='Send message to admin.')
    p_syn.add_argument('--send-msg-to-user', action='store_true',
                       help='Send message to user.')
    p_syn.add_argument('--period', default=1800,
                       help='Time period in second to download submissions from synapse '
                            'and score them')
    p_syn.add_argument('--admin-id', nargs='+', default=['3345120'],
                       help='Admin\'s Synapse ID (as string) ')
    args = parser.parse_args()

    return args

def validate_submission(submission, status, args, syn):
    status['status'] = 'INVALID'

    submission_dir = os.path.join(
        os.path.abspath(args.submission_dir), submission.id)
    mkdir_p(submission_dir)

    metadata = {
        'id': submission.id,
        'team': 'undefined'
    }

    try:
        metadata['team'] = get_team_name(syn, None, submission.teamId)

        log.info('Downloading submission... {}'.format(submission.id))
        submission = syn.getSubmission(
            submission, 
            downloadLocation=submission_dir, 
            ifcollision='overwrite.local'
        )
        print()
        submission_fname = submission.filePath
        cell, assay = parse_submission_filename(submission_fname)
        if not is_valid_round2_cell_assay(cell, assay):
            raise Exception('Invalid cell/assay combination for '
                            'round2 round')

        log.info('Downloading done {}, {}, {}, {}, {}'.format(
            submission_fname, submission.id,
            submission.teamId, cell, assay))

        # read pred npy (submission)
        log.info('Validating bigwig...{}'.format(submission.id))
        valid, message = validate(submission_fname, args.window_size)

        if valid:
            status['status'] = 'VALIDATED'
            subject = 'Successfully validated submission %s %s %s:\n' % (
                submission.name, submission.id, submission.teamId)
        else:
            subject = 'Invalid submission %s %s %s:\n' % (
                submission.name, submission.id, submission.teamId)
            # delete file
            shutil.rmtree(submission_dir)

        log.info(subject + message)

    except Exception as ex1:
        if 'teamId' in submission:
            teamId = submission.teamId
        else:
            teamId = 'undefined'
        subject = 'Error validating submission %s %s %s:\n' % (
            submission.name, submission.id, teamId)
        st = StringIO()
        traceback.print_exc(file=st)
        message = st.getvalue()
        log.error(subject + message)
        # delete file
        shutil.rmtree(submission_dir)

    # send message
    users_to_send_msg = []
    if args.send_msg_to_user:
        users_to_send_msg.append(submission.userId)
    if args.send_msg_to_admin:
        users_to_send_msg.extend(args.admin_id)
    send_message(syn, users_to_send_msg, subject, message)

    if not args.dry_run:
        status['annotations'] = synapseclient.annotations.to_submission_status_annotations(
            metadata, is_private=False)
        status = syn.store(status)

    return status


def main():
    log.info('Parsing arguments...')

    args = parse_arguments()

    if args.team_name_tsv is not None:
        team_name_dict = parse_team_name_tsv(args.team_name_tsv)
    else:
        team_name_dict = None
    print(team_name_dict)

    syn = synapseclient.login()
    t0 = time.perf_counter()

    while True:
        try:
            if not args.update_wiki_only:
                evaluation = syn.getEvaluation(args.eval_queue_id)

                # init multiprocessing
                pool = multiprocessing.Pool(args.nth)

                # distribute jobs
                ret_vals = []
                for submission, status in syn.getSubmissionBundles(evaluation, status='RECEIVED'):
                    ret_vals.append(
                        pool.apply_async(validate_submission,
                                         (submission, status, args, syn,
                                          gene_annotations, enh_annotations)))
                # gather
                for r in ret_vals:
                    r.get(BIG_INT)

                pool.close()
                pool.join()

            update_wiki(syn, team_name_dict, args)

        except Exception as ex1:
            st = StringIO()
            traceback.print_exc(file=st)
            message = st.getvalue()

            subject = 'Server error:'
            if args.send_msg_to_admin:
                send_message(syn, args.admin_id, subject, message)
            log.error(message)

        log.info('Waiting for new submissions...')
        while time.perf_counter() - t0 < args.period:
            time.sleep(60)
        t0 = time.perf_counter()

    log.info('All done')
 

if __name__ == '__main__':
    main()

