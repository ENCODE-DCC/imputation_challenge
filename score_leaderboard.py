#!/usr/bin/env python3
"""Imputation challenge scoring script for leaderboard
Author:
    Jin Lee (leepc12@gmail.com)
"""

import os
import re
import time
import shutil
#import gc
import traceback
import synapseclient
import multiprocessing
from bw_to_npy import load_bed, load_npy, bw_to_dict
from score import parse_submission_filename, score
from score_metrics import Score
from rank import calc_global_ranks, get_cell_name, get_assay_name, get_team_name, parse_team_name_tsv
from db import write_to_db, ScoreDBRecord, DB_QUERY_GET, read_scores_from_db
from io import StringIO
from logger import log


BIG_INT = 99999999  # for multiprocessing


LEADERBOARD_ROUND_VALID_CELL_ASSAY = [
'C02M22',
'C03M02',
'C04M16',
'C09M20',
'C10M17',
'C12M16',
'C12M32',
'C13M20',
'C16M17',
'C17M04',
'C17M19',
'C17M29',
'C17M32',
'C18M21',
'C18M25',
'C20M22',
'C23M03',
'C23M07',
'C23M26',
'C23M34',
'C24M17',
'C24M25',
'C25M21',
'C25M26',
'C27M03',
'C27M13',
'C27M24',
'C27M26',
'C29M29',
'C31M25',
'C32M08',
'C32M12',
'C32M20',
'C34M02',
'C34M32',
'C36M18',
'C37M29',
'C45M22',
'C46M10',
'C46M18',
'C46M21',
'C46M35',
'C47M18',
'C48M16',
'C50M02'
]

def is_valid_leaderboard_cell_assay(cell, assay):
    return (cell + assay) in LEADERBOARD_ROUND_VALID_CELL_ASSAY


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


# markdown supertable to show submission status
WIKI_TEMPLATE_SUBMISSION_STATUS = \
'${{supertable?path=%2Fevaluation%2Fsubmission%2Fquery%3Fquery%3D\
select+%2A+from+evaluation_{eval_queue_id}+\
&paging=true&queryTableResults=true&showIfLoggedInOnly=false&pageSize=100\
&showRowNumber=false&jsonResultsKeyName=rows\
&columnConfig0=none%2CID%2CobjectId%3B%2CNONE\
&columnConfig1=none%2CteamId%2CteamId%3B%2CNONE\
&columnConfig2=none%2Cteam%2Cteam%3B%2CNONE\
&columnConfig3=epochdate%2CDate%2CcreatedOn%3B%2CNONE\
&columnConfig4=none%2Cname%2Cname%3B%2CNONE\
&columnConfig5=none%2Cstatus%2Cstatus%3B%2CNONE%2C4\
}}\n\n'

# markdown supertable to show submission status
WIKI_TEMPLATE_SUBMISSION_SCORE = \
'${{supertable?path=%2Fevaluation%2Fsubmission%2Fquery%3Fquery%3D\
select+%2A+from+evaluation_{eval_queue_id}+\
where%2Bstatus%3D%3D%2522SCORED%2522%2Band%2Bcell%3D%3D%2522{cell}%2522%2B\
and%2Bassay%3D%3D%2522{assay}%2522%2B\
&paging=true&queryTableResults=true&showIfLoggedInOnly=false&pageSize=100\
&showRowNumber=false&jsonResultsKeyName=rows\
&columnConfig0=none%2CID%2CobjectId%3B%2CNONE\
&columnConfig1=none%2CteamId%2CteamId%3B%2CNONE\
&columnConfig2=none%2Cteam%2Cteam%3B%2CNONE\
&columnConfig3=epochdate%2CDate%2CcreatedOn%3B%2CNONE\
&columnConfig4=none%2Cmse%2Cmse%3B%2CNONE\
&columnConfig5=none%2Cgwcorr%2Cgwcorr%3B%2CNONE\
&columnConfig6=none%2Cgwspear%2Cgwspear%3B%2CNONE\
&columnConfig7=none%2Cmseprom%2Cmseprom%3B%2CNONE\
&columnConfig8=none%2Cmsegene%2Cmsegene%3B%2CNONE\
&columnConfig9=none%2Cmseenh%2Cmseenh%3B%2CNONE\
&columnConfig10=none%2Cmsevar%2Cmsevar%3B%2CNONE\
&columnConfig11=none%2Cmse1obs%2Cmse1obs%3B%2CNONE\
&columnConfig12=none%2Cmse1imp%2Cmse1imp%3B%2CNONE\
}}\n\n'

RE_PATTERN_SUBMISSION_FNAME = r'^C\d\dM\d\d.*(bw|bigwig|bigWig|BigWig)'

def update_wiki(syn, team_name_dict, args):
    # calculate ranks and update leaderboard wiki
    log.info('Updating wiki...')
    rows = read_scores_from_db(args.db_file, args.chrom)
    _, _, markdown_per_cell_assay, markdown_overall = calc_global_ranks(
        rows, args.measures_to_use, team_name_dict, syn)

    wiki_id_map = {
        k.split(':')[0]: k.split(':')[1] for k in args.leaderboard_wiki_id.split(',')
    }
    for k, wiki_id in wiki_id_map.items():
        w = syn.getWiki(args.project_id, wiki_id)

        if k == 'submission_status':
            title = 'Submission status'
            markdown = WIKI_TEMPLATE_SUBMISSION_STATUS.format(
                eval_queue_id=args.eval_queue_id)

        elif k == 'overall':
            title = 'Overall ranks'
            markdown = markdown_overall

        elif k.startswith('C'):
            title = '{} {}'.format(k, get_cell_name(k))
            markdown = ''
            if k in markdown_per_cell_assay:
                for assay, m in markdown_per_cell_assay[k].items():
                    markdown += m + '\n'
                    markdown += WIKI_TEMPLATE_SUBMISSION_SCORE.format(
                        eval_queue_id=args.eval_queue_id,
                        cell=k, assay=assay)
        else:
            continue

        w.markdown = markdown
        w.title = title
        w = syn.store(w)

    return None

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
    p_sys.add_argument('--team-name-tsv',
                        help='TSV file with team_id/team_name (1st col/2nd col).')
    p_syn = parser.add_argument_group(
                        title='Communitation with synapse')
    p_syn.add_argument('--dry-run', action='store_true',
                       help='Do not update submission\'s status on synapse.')
    p_syn.add_argument('--project-id', default='syn17083203',
                       help='Synapse project ID.')
    p_syn.add_argument('--leaderboard-wiki-id',
                       default='overall:594046,submission_status:594047,'
                               'C02:594048,C03:594049,C04:594050,C09:594051,C10:594052,'
                               'C12:594053,C13:594054,C16:594055,C17:594056,C18:594057,'
                               'C20:594058,C23:594059,C24:594060,C25:594061,C27:594062,'
                               'C29:594063,C31:594064,C32:594065,C34:594068,C36:594069,'
                               'C37:594070,C45:594071,C46:594072,C47:594073,C48:594074,'
                               'C50:594075',
                       help='Comma-delimited Synapse wiki ID for leaderboard. '
                            'Required items: overall, submission_status, C??'
                            'Format example: "overall:594046,'
                            'submission_status:594047"')
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
    status['status'] = 'INVALID'

    submission_dir = os.path.join(
        os.path.abspath(args.submission_dir), submission.id)
    mkdir_p(submission_dir)

    chosen_score = None  # first bootstrap score
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
        if not is_valid_leaderboard_cell_assay(cell, assay):
            raise Exception('Invalid cell/assay combination for '
                            'leaderboard round')

        log.info('Downloading done {}, {}, {}, {}, {}'.format(
            submission_fname, submission.id,
            submission.teamId, cell, assay))

        # read pred npy (submission)
        log.info('Converting to dict...{}'.format(submission.id))
        y_pred_dict = bw_to_dict(submission_fname, args.chrom,
                                 args.window_size, args.blacklist_file,
                                 args.validated)
        #gc.collect()
        # read truth npy
        npy_true = os.path.join(
            args.true_npy_dir,
            '{}{}.npy'.format(cell, assay))
        y_true_dict = bw_to_dict(npy_true, args.chrom,
                                 args.window_size, args.blacklist_file)
        #gc.collect()
        # read var npy
        if args.var_npy_dir is not None:   
            var_npy = os.path.join(
                args.var_npy_dir,
                'var_{}.npy'.format(assay))
            y_var_dict = load_npy(var_npy)
        else:
            y_var_dict = None
        #gc.collect()

        score_outputs = []
        for k, bootstrap_chrom in args.bootstrap_chrom:
            # score it for each bootstrap chroms
            log.info('Scoring... k={}, submission_id={}'.format(k, submission.id))
            r = score(y_pred_dict, y_true_dict, bootstrap_chrom,
                      gene_annotations, enh_annotations,
                      args.window_size, args.prom_loc,
                      y_var_dict)
            #gc.collect()  # free memory for bootstrapped arrays
            log.info('Scored: {}, {}, {}, {}'.format(
                submission.id, submission.teamId, k, r))
            score_outputs.append((k, r))

        # score to be shown on wiki (first bootstrap score)
        chosen_score = score_outputs[0][1]

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
        status['status'] = 'SCORED'

        # free memory
        y_pred_dict = None
        y_true_dict = None
        y_var_dict = None
        #gc.collect()

        subject = 'Successfully scored submission %s %s %s:\n' % (
            submission.name, submission.id, submission.teamId)
        message = 'Score (bootstrap_idx: score)\n'
        message += '\n'.join(
            ['{}: {}'.format(k, s) for k, s in score_outputs])
        log.info(subject + message)

    except Exception as ex1:
        if 'teamId' in submission:
            teamId = submission.teamId
        else:
            teamId = 'undefined'
        subject = 'Error scoring submission %s %s %s:\n' % (
            submission.name, submission.id, teamId)
        st = StringIO()
        traceback.print_exc(file=st)
        message = st.getvalue()
        log.error(subject + message)

    finally:
        # remove submissions (both bigwig, npy) to save disk space
        shutil.rmtree(submission_dir)
        pass

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

    enh_annotations = load_bed(args.enh_annotations)
    gene_annotations = load_bed(args.gene_annotations)

    # do GC manually
    #gc.disable()

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
                        pool.apply_async(score_submission,
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

