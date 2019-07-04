#!/usr/bin/env python3
"""Imputation challenge synapse interface

Author:
    Jin Lee (leepc12@gmail.com)
"""

import synapseclient
import exceptions
from logger import log

syn = None

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


def init_synapse():
    """Login on synapse. Authentication must be done in CLI
    """
    global syn
    syn = synapseclient.login()    


def get_team_name(team_id):
    global syn

    if syn is None:
        return TEAM_ID_FOR_ROUND1[team_id]
    else:
        team = syn.restGET('/team/{id}'.format(id=team_id))
        return team['name']


def download_submissions(eval_queue_id):
    """Download submission from an evaluation queue
    """
    global syn
    result = {}

    assert(syn is not None)
    evaluation = syn.getEvaluation(eval_queue_id)
    for submission, status in syn.getSubmissionBundles(
        evaluation, status='RECEIVED'):

        try:
            ## refetch the submission so that we get the file path
            ## to be later replaced by a "downloadFiles" flag on getSubmissionBundles
            submission_id = submission.id
            team_id = submission.teamId
            
            print(submission_id, team_id, submission)
            submission = syn.getSubmission(submission, downloadLocation='tmp/{}'.format(submission_id))
            # print(submission.entityId, submission.name, submission.evaluationId)
        except:
            mark_submission_status(submission_id, 'FAILED_TO_DL')

    return result


def mark_submission_status(submission_id, status):
    