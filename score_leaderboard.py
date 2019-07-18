#!/usr/bin/env python3
"""Imputation challenge scoring script for leaderboard
Author:
    Jin Lee (leepc12@gmail.com)
"""

import os
import time
import shutil
import gc
import traceback
import synapseclient
import multiprocessing
from bw_to_npy import load_bed, load_npy, bw_to_dict
from score import parse_submission_filename, score
from rank import calc_global_ranks, get_cell_name, get_assay_name, get_team_name
from db import write_to_db, ScoreDBRecord, DB_QUERY_GET, read_scores_from_db
from io import StringIO
from logger import log


BIG_INT = 99999999  # for multiprocessing


def mkdir_p(path):
    import errno    
    import os

    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def send_message(syn, user_ids, subject, message):
    if len(user_ids) == 0:
        return None
    try:
        response = syn.sendMessage(
            userIds=user_ids,
            messageSubject=subject,
            messageBody=message,
            contentType="text/html")
        log.info("Message sent: {}".format(str(response).encode('utf-8')))
    except Exception as ex0:
        log.error("Failed to send message: {}".format(ex0))
        response = None

    return response


def update_wiki(syn, project_id, wiki_id, markdown):
    w = syn.getWiki(project_id, wiki_id)
    w.markdown = markdown
    w = syn.store(w)
    return w


def parse_arguments():
    import argparse
    import os

    py_path = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser(
        description='ENCODE Imputation Challenge leaderboard script. ')
    parser.add_argument('eval_queue_id',
                        help='Synapse evaluation queue ID to retreive submissions from.')
    parser.add_argument('true_npy_dir',
                        help='Directory for truth .npy files. '
                             'All .npy files will be used for scoring against the submission')
    parser.add_argument('--var-npy-dir',
                        help='Directory for var .npy files. '
                             'All var_CXX.npy files will be used for scoring.')
    p_score = parser.add_argument_group(
                        title='Scoring parameters')
    p_score.add_argument('--chrom', nargs='+',
                         default=['all'],
                         help='List of chromosomes to be combined to be '
                              'scored. '
                              'Set as "all" (default) to score for all '
                              'chromosomes. '
                              '(e.g. "all" or "chr3 chr21") '
                              'It should be "all" to write scores to DB file')
    p_score.add_argument('--bootstrap-chrom', nargs='*', default=[],
                         help='Bootstrapped chromosome groups. '
                              'Delimiter is whitespace for groups and '
                              'comma(,) in each group. Order is important.'
                              'e.g. "chr1,chr2 chr1,chrX chr2,chrX" means '
                              'three groups: (chr1,chr2), (chr1,chrX), (chr2,chrX)')
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
    p_score.add_argument('--blacklist-file',
                         default=os.path.join(
                            py_path,
                            'annot/hg38/hg38.blacklist.bed.gz'),
                         help='Blacklist BED file. Bootstrap label will be '
                              'generated after removing overlapping regions '
                              'defined in this file.')
    p_score.add_argument('--window-size', default=25, type=int,
                         help='Window size for bigwig in bp')
    p_score.add_argument('--prom-loc', default=80, type=int,
                         help='Promoter location in a unit of window size '
                              '(--window-size). This is not in bp')
    p_score.add_argument('--measures-to-use', nargs='+',
                        default=['mse', 'gwcorr', 'gwspear', 'mseprom',
                                 'msegene', 'mseenh', 'msevar', 'mse1obs',
                                 'mse1imp'],
                        help='List of performance measures to be used for ranking')
    p_score.add_argument('--validated', action='store_true',
                         help='For validated submissions '
                              'with fixed interval length of 25 and valid '
                              'chromosome lengths. It will skip interpolation')
    p_score.add_argument('--update-wiki-only', action='store_true',
                         help='Update wiki based on DB file (--db-file) without '
                              'scoring submissions')
    #p_score.add_argument('--normalize-with-robust-min-max', action='store_true',
    #                     help='Normalize with robust min max.')
    p_out = parser.add_argument_group(
                        title='Output database file')
    p_out.add_argument('--db-file',
                       help='Write metadata/scores to SQLite DB file')
    p_sys = parser.add_argument_group(
                        title='System and resource settings')
    p_sys.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize scoring (per) submission')
    p_syn = parser.add_argument_group(
                        title='Communitation with synapse')
    p_syn.add_argument('--dry-run', action='store_true',
                       help='Do not update submission\'s status on synapse.')
    p_syn.add_argument('--project-id', default='syn17083203',
                       help='Synapse project ID.')
    p_syn.add_argument('--leaderboard-wiki-id', default='parent:594012',
                       help='Comma-delimited Synapse wiki ID for leaderboard. '
                            'Format example: "parent:594012,C01:594013,C02:594014"')
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

    if args.db_file is not None:
        args.db_file = args.db_file.strip("'")
        if not os.path.exists(args.db_file):
            raise ValueError('DB file does not exists')

    if args.chrom == ['all']:
        args.chrom = ['chr' + str(i) for i in range(1, 23)] + ['chrX']
    args.chrom = sorted(args.chrom)

    if len(args.bootstrap_chrom) == 0:
        args.bootstrap_chrom = [(-1, args.chrom)]
    else:
        for i, _ in enumerate(args.bootstrap_chrom):
            args.bootstrap_chrom[i] = (i, args.bootstrap_chrom[i].split(','))
    log.info(args.bootstrap_chrom)

    return args


def score_submission(submission, status, args, syn,
                     gene_annotations, enh_annotations):
    status.status = "INVALID"
    submission_dir = os.path.join(
        os.path.abspath(args.submission_dir), submission.id)
    mkdir_p(submission_dir)

    chosen_score = None  # first bootstrap score
    metadata = {
        'id': submission.id,
        'team': get_team_name(syn, None, submission.teamId)
    }

    try:
        log.info('Downloading submission... {}'.format(submission.id))
        submission = syn.getSubmission(
            submission, 
            downloadLocation=submission_dir, 
            ifcollision='overwrite.local'
        )
        print()
        submission_fname = submission.filePath
        cell, assay = parse_submission_filename(submission_fname)
        log.info('Downloading done {}, {}, {}, {}, {}'.format(
            submission_fname, submission.id,
            submission.teamId, cell, assay))

        # read pred npy (submission)
        log.info('Converting to dict...{}'.format(submission.id))
        y_pred_dict = bw_to_dict(submission_fname, args.chrom,
                                 args.window_size, args.blacklist_file,
                                 args.validated)
        gc.collect()
        # read truth npy
        npy_true = os.path.join(
            args.true_npy_dir,
            '{}{}.npy'.format(cell, assay))
        y_true_dict = bw_to_dict(npy_true, args.chrom,
                                 args.window_size, args.blacklist_file)
        gc.collect()
        # read var npy
        if args.var_npy_dir is not None:   
            var_npy = os.path.join(
                args.var_npy_dir,
                'var_{}.npy'.format(assay))
            y_var_dict = load_npy(var_npy)
        else:
            y_var_dict = None
        gc.collect()

        score_outputs = []
        for k, bootstrap_chrom in args.bootstrap_chrom:
            # score it for each bootstrap chroms
            log.info('Scoring... k={}, submission_id={}'.format(k, submission.id))
            r = score(y_pred_dict, y_true_dict, bootstrap_chrom,
                      gene_annotations, enh_annotations,
                      args.window_size, args.prom_loc,
                      y_var_dict)
            gc.collect()  # free memory for bootstrapped arrays
            log.info('Scored: {}, {}, {}, {}'.format(
                submission.id, submission.teamId, k, r))
            score_outputs.append((k, r))

        chosen_score = score_outputs[0]

        # write to db and report
        for k, score_output in score_outputs:
            if not args.dry_run:
                score_db_record = ScoreDBRecord(
                    int(submission.id),
                    int(submission.teamId),
                    submission_fname,
                    cell,
                    assay,
                    k,
                    *score_output)
                write_to_db(score_db_record, args.db_file)
        # mark is as scored
        status.status = "SCORED"

        # free memory
        y_pred_dict = None
        y_true_dict = None
        y_var_dict = None
        gc.collect()

        subject = 'Successfully scored submission %s %s %s:\n' % (
            submission.name, submission.id, submission.userId)
        message = 'Score (bootstrap_idx: score)\n'
        message += '\n'.join(
            ['{}: {}'.format(k, s) for k, s in score_outputs])
        log.info(subject + message)

    except Exception as ex1:
        subject = 'Error scoring submission %s %s %s:\n' % (
            submission.name, submission.id, submission.userId)
        st = StringIO()
        traceback.print_exc(file=st)
        message = st.getvalue()
        log.error(subject + message)

    finally:
        # remove submissions (both bigwig, npy) to save disk space
        # shutil.rmtree(submission_dir)

    # send message
    users_to_send_msg = []
    if args.send_msg_to_user:
        users_to_send_msg.append(submission.userId)
    if args.send_msg_to_admin:
        users_to_send_msg.extend(args.admin_id)
    send_message(syn, users_to_send_msg, subject, message)

    if not args.dry_run:
        # update metadata with score
        if chosen_score is not None:
            for field in Score._fields:
                metadata[field] = getattr(chosen_score, field)
            metadata['cell'] = cell
            metadata['assay'] = assay

        status.annotations = synapseclient.annotations.to_submission_status_annotations(
            metadata, is_private=False)
        status = syn.store(status)

    return None


def main():
    log.info('Parsing arguments...')

    args = parse_arguments()

    enh_annotations = load_bed(args.enh_annotations)
    gene_annotations = load_bed(args.gene_annotations)

    # do GC manually
    gc.disable()

    syn = synapseclient.login()

    while True:
        try:
            if not args.update_wiki_only:
                evaluation = syn.getEvaluation(args.eval_queue_id)

                # init multiprocessing
                pool = multiprocessing.Pool(args.nth)

                # distribute jobs
                ret_vals = []
                #for submission, status in syn.getSubmissionBundles(evaluation, status='RECEIVED'):
                for submission, status in syn.getSubmissionBundles(evaluation):
                    #print(submission, status)
                    #if status == 'SCORED':
                    #    continue
                    ret_vals.append(
                        pool.apply_async(score_submission,
                                         (submission, status, args, syn,
                                          gene_annotations, enh_annotations)))
                # gather
                for r in ret_vals:
                    r.get(BIG_INT)

                pool.close()
                pool.join()

            # calculate ranks and update leaderboard wiki
            log.info('Updating ranks...')
            rows = read_scores_from_db(args.db_file, args.chrom)
            _, _, markdown = calc_global_ranks(
                rows, args.measures_to_use, None, syn)

            update_wiki(syn, args.project_id, args.leaderboard_wiki_id,
                        markdown)

        except Exception as ex1:
            st = StringIO()
            traceback.print_exc(file=st)
            message = st.getvalue()

            subject = 'Server error:'
            if args.send_msg_to_admin:
                send_message(syn, args.admin_id, subject, message)
            log.error(message)

        log.info('Waiting for {} secs to check new submissions on '
                 'syn eval queue'.format(args.period))
        time.sleep(args.period)

    log.info('All done')
 

if __name__ == '__main__':
    main()

