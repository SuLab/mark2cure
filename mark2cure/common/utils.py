from django.conf import settings
from django.db.models import Count
from django.core.mail import send_mail

from mark2cure.document.models import Activity

from random import shuffle
import logging
logger = logging.getLogger(__name__)


def experiment_gm_routing(user, doc_arr):
    #user_profile = user.userprofile
    '''
      I need to figure out which of the current golden
      documents have been done the least
      1. Array of all document_id options
      2. Array of all documents already done by that user
      3. Sorting by completion counts
    '''
    #prev_docs = Activity.objects.filter(user = user, experiment = settings.EXPERIMENT if user.userprofile.mturk else None).values_list('document__pk', flat = True).all()
    shuffle(doc_arr)
    return doc_arr[0]


def experiment_routing(user, doc_arr, n_value = 15):
    user_profile = user.userprofile
    '''
      I need to figure out which of the current experiment
      documents are still available for work to be done on
      1. Array of all document_id options
      2. Array of all documents already done by that user
      3. Array of all documents from options which have already been completed their max times
      4. Array of (1 - 2) - 3
    '''

    prev_docs = Activity.objects.filter(user = user, experiment = settings.EXPERIMENT if user.userprofile.mturk else None).exclude(user__userprofile__ignore = True).values_list('document__pk', flat = True).all()

    # This is actually very fast, returns 593 doc ids
    #experiment_docs = Document.objects.filter(source = 'NCBI_corpus_training').values_list('pk', flat = True).all()
    # Make a copy so we can remove from copy
    experiment_docs = list(doc_arr)
    for x in prev_docs:
        if x in experiment_docs: experiment_docs.remove(x)

    '''
      Count the number of times the experimental docs have been completed for Concept Recognition
      within this Experiment (if relevant), that aren't GM (none should be anyway)
    '''
    activities = Activity.objects.filter(
        document__pk__in = experiment_docs,
        task_type = 'cr',
        experiment = settings.EXPERIMENT if user.userprofile.mturk else None).exclude(submission_type = 'gm', user__userprofile__ignore = True).values('document').annotate(Count('document'))

    experiment_docs_completed = [item['document'] for item in activities if item['document__count'] >= n_value]

    # Email us to let us know when the K saturates on these
    #if len(experiment_docs_completed) > 15:
    #    send_mail('[Mark2Cure] Document Completion Milestone #{0}'.format(settings.EXPERIMENT),
    #            '{0} documents have had all their work completed.'.format(len(experiment_docs_completed)),
    #            settings.SERVER_EMAIL,
    #            [email[1] for email in settings.MANAGERS])

    for x in experiment_docs_completed:
        if x in experiment_docs: experiment_docs.remove(x)

    logger.debug("Experiment Docs remaining count {0} for {1} in experiment {2}".format(len(experiment_docs), user.username, settings.EXPERIMENT))

    shuffle(experiment_docs)

    if len(experiment_docs) == 0:
        # All the docs have been completed (allow unlimited HITs with extra protection)
        return False
    else:
        return experiment_docs[0]


