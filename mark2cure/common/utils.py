from django.conf import settings
from django.db.models import Count

from mark2cure.document.models import Activity

from datetime import datetime
import datetime, random


def experiment_routing(user, n_count, gm_occurance = 4, k_max = 3):
    user_profile = user.userprofile

    # (TODO) Pretect from out of range & smarter gold injection
    gm_docs = [282, 310, 363, 326, 322]
    if n_count % gm_occurance == 0:
        gm_index = (n_count / gm_occurance) - 1
        if (gm_index+1) > len(gm_docs):
            return gm_docs[gm_index]

    if user_profile.mturk:
        prev_docs = Activity.objects.filter(user=user, experiment=settings.EXPERIMENT).values('document__pk', flat=True).all()
    else:
        prev_docs = Activity.objects.filter(user=user).values_list('document__pk', flat=True).all()

    '''
      I need to figure out which of the current experiment
      documents are still available for work to be done on
      1. Array of all document_id options
      2. Array of all documents already done by that user
      3. Array of all documents from options which have already been completed their max times
      4. Array of (1 - 2) - 3
    '''
    experiment_docs = [2787, 3195, 3357, 2726, 2637, 3030, 3203, 3314, 3077, 2369, 2394, 3003, 3567, 3166, 3177, 2152, 2661, 2236, 2193, 2878]
    for x in prev_docs:
      experiment_docs.remove(x)

    activities = Activity.objects.filter(document__pk__in=experiment_docs, task_type='cr', submission_type='gm', experiment=None).values('document').annotate(Count('document'))
    experiment_docs_completed = [item['document'] for item in activities if item['document__count'] >= k_max]

    for x in experiment_docs_completed:
        experiment_docs.remove(x)


    random.shuffle(experiment_docs)
    return experiment_docs[0]




