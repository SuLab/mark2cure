from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import Count
from django.core.mail import send_mail

from mark2cure.document.models import Activity

from datetime import datetime
import datetime, random, logging
logger = logging.getLogger(__name__)


def experiment_routing(user, n_count, k_max = 8):
    user_profile = user.userprofile

    gm_dict = {
        9: 282,
        17: 310,
        40: 363,
        41: 326,
        42: 322
        }

    if n_count in gm_dict:
        return gm_dict[n_count]

    '''
      I need to figure out which of the current experiment
      documents are still available for work to be done on
      1. Array of all document_id options
      2. Array of all documents already done by that user
      3. Array of all documents from options which have already been completed their max times
      4. Array of (1 - 2) - 3
    '''
    prev_docs = Activity.objects.filter(user=user, experiment= settings.EXPERIMENT if user.userprofile.mturk else None).values_list('document__pk', flat=True).all()
    # Make a copy so we can remove from copy
    experiment_docs = settings.EXPERIMENT_DOCS
    for x in prev_docs:
        if x in experiment_docs: experiment_docs.remove(x)

    '''
      Count the number of times the experimental docs have been completed for Concept Recognition
      within this Experiment (if relevant), that aren't GM (none should be anyway)
    '''
    activities = Activity.objects.filter(
        document__pk__in=experiment_docs,
        task_type='cr',
        experiment= settings.EXPERIMENT if user.userprofile.mturk else None).exclude(submission_type='gm').values('document').annotate(Count('document'))
    experiment_docs_completed = [item['document'] for item in activities if item['document__count'] >= k_max]

    for x in experiment_docs_completed:
        if x in experiment_docs: experiment_docs.remove(x)

    logger.debug("Experiment Docs remaining count {0} for {1} in experiment {2}".format(len(experiment_docs), user.username, settings.EXPERIMENT))

    random.shuffle(experiment_docs)

    if len(experiment_docs) < 5:
        send_mail('[Mark2Cure] Experiment Doc Remaining #{0}'.format(settings.EXPERIMENT),
                'There are only {0} documents that have work left to be done.'.format(len(experiment_docs)),
                settings.SERVER_EMAIL,
                [email[1] for email in settings.MANAGERS])

    if len(experiment_docs) == 0:
        # All the docs have been completed (allow unlimited HITs with extra protection)
        return False
    else:
        return experiment_docs[0]


