#!/usr/bin/env python3
"""Imputation challenge synapse test script
Author:
    Jin Lee (leepc12@gmail.com)
"""

import synapseclient
# from challenge_metadata import init_synapse, download_submissions_from_syn_eval_queue

eval_queue_id = 9614263

syn = synapseclient.login()
evaluation = syn.getEvaluation(eval_queue_id)
print(evaluation)

for submission, status in syn.getSubmissionBundles(
        evaluation, status='RECEIVED'):

    ## refetch the submission so that we get the file path
    ## to be later replaced by a "downloadFiles" flag on getSubmissionBundles
    submission_id = submission.id
    team_id = submission.teamId
    
    print(submission_id, team_id, submission)
    submission = syn.getSubmission(submission, downloadLocation='tmp/{}'.format(submission_id))
    # print(submission.entityId, submission.name, submission.evaluationId)

    # print "validating", submission.id, submission.name
    # try:
    #     is_valid, validation_message = conf.validate_submission(
    #         evaluation, submission, ifcollision='overwrite.local')
    # except Exception as ex1:
    #     is_valid = False
    #     print "Exception during validation:", type(ex1), ex1, ex1.message
    #     traceback.print_exc()
    #     validation_message = str(ex1)

    # status.status = "VALIDATED" if is_valid else "INVALID"
