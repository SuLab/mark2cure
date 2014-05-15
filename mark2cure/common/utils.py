from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import Count

from mark2cure.document.models import Activity

from datetime import datetime
import datetime, random


def experiment_routing(user, n_count, gm_occurance = 4, k_max = 3):
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

    # if n_count % gm_occurance == 0:
    #     gm_index = (n_count / gm_occurance) - 1
    #     if (gm_index+1) > len(gm_docs):
    #         return gm_docs[gm_index]

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

    activities = Activity.objects.filter(document__pk__in=experiment_docs, task_type='cr', submission_type='gm', experiment=None).values('document').annotate(Count('document'))
    experiment_docs_completed = [item['document'] for item in activities if item['document__count'] >= k_max]

    for x in experiment_docs_completed:
        experiment_docs.remove(x)


    random.shuffle(experiment_docs)
    return experiment_docs[0]




