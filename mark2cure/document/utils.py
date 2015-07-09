from django.conf import settings
from django.core.exceptions import ObjectDoesNotExist

from .models import Document, Section, View, Annotation

from django.contrib.auth.models import User

from datetime import datetime
from Bio import Entrez, Medline

import random


def select_best_opponent(task, document, player):
    '''
        Select the best opponate for a certain task meaning a
        certain document scoped to a quest for a certain task type.

        1) First weight by GM, if one is available always prefer it over
            other users

        2) Select users with non-empty response for this View (Document scoped to Quest)
            * Empty checks for len() > 0

            //-- Stop here no need to go into this complexity for Experiment 2
            2.1) Sort by internal F score average, select top 3 (over threshold), pick 1 at random
                * If score is not present, calcuate it for each user to be cached for next time
                * Score is average F score on last 10 computable documents
                  * This creates a rolling window, could do # or time based.
                    * Bad if documents completed long time ago for # based

            2.2) Save new score for current comparision
                * Save score if compared to a Golden Master document (current pool size is limited)
                * Save score if compared to a top ranked (1 of the 3) that was randomly selected
                * Do not have novel annotations or comparisions to high skill badged people

        Downfalls:
            * Good user could annotate all over the place, generate
              random noise but would pass validation as they're top ranked and len() > 0
            * Does not account for total game play (this might be good), uses rolling window
              of history
    '''
    # Select all (**including uncompleted**) other user started quests
    # that have been completed
    others_quest_relationships = task.userquestrelationship_set.exclude(user=player)

    # If the known GM User is in the DB, use them for partner comparison
    gm_user_query = User.objects.filter(username='GATTACA')
    gm_user = None
    if gm_user_query.exists():
        gm_user = gm_user_query.first()
        if others_quest_relationships.exists() and \
                others_quest_relationships.filter(user=gm_user).exists() and \
                others_quest_relationships.filter(user=gm_user).first().views.filter(section__document=document, completed=True).exists():
            # There is an "expert's" annotations (GM) so
            # show those as the partner's
            return gm_user

    # Gather users from completed documents
    # that may come from uncompleted quests
    previous_users = []
    for quest_relationship in others_quest_relationships.exclude(user=gm_user, user__pk__in=[107,]):
        if quest_relationship.views.filter(section__document=document, completed=True).exists():
            # (TODO) Don't add option of them unless they've submitted 1+ Annotations
            view_ids = quest_relationship.views.filter(section__document=document, completed=True).values_list('pk', flat=True)
            if Annotation.objects.filter(view__pk__in=view_ids).exists():
                previous_users.append(quest_relationship.user)

    if others_quest_relationships.exists() and len(previous_users):
        # No expert around so select a previous user at random
        random.shuffle(previous_users)
        selected_user = previous_users[0]
        return selected_user

    return None


def determine_f(true_positive, false_positive, false_negative):
    if float(true_positive + false_positive) == 0.0:
        return (0.0, 0.0, 0.0)

    if float(true_positive + false_negative) == 0.0:
        return (0.0, 0.0, 0.0)

    precision = true_positive / float(true_positive + false_positive)
    recall = true_positive / float(true_positive + false_negative)

    if float(precision + recall) > 0.0:
        f = (2 * precision * recall) / (precision + recall)
        return (precision, recall, f)
    else:
        return (0.0, 0.0, 0.0)


def match_exact(gm_ann, user_anns):
    for user_ann in user_anns:
        if user_ann.is_exact_match(gm_ann):
            return True
    return False


def generate_results(user_views, gm_views):
    '''
      This calculates the comparsion overlap between two arrays of dictionary terms

      It considers both the precision p and the recall r of the test to compute the score:
      p is the number of correct results divided by the number of all returned results
      r is the number of correct results divided by the number of results that should have been returned.
      The F1 score can be interpreted as a weighted average of the precision and recall, where an F1 score reaches its best value at 1 and worst score at 0.

     tp  fp
     fn  *tn

    '''
    gm_annotations = Annotation.objects.filter(view__pk__in=[v.pk for v in gm_views if type(v) is View])
    user_annotations = Annotation.objects.filter(view__pk__in=[v.pk for v in user_views if type(v) is View])

    true_positives = [gm_ann for gm_ann in gm_annotations if match_exact(gm_ann, user_annotations)]

    # Annotations the user submitted that were wrong (the User set without their True Positives)
    # false_positives = user_annotations - true_positives
    false_positives = user_annotations
    for tp in true_positives:
        false_positives = false_positives.exclude(start=tp.start, text=tp.text)

    # # Annotations the user missed (the GM set without their True Positives)
    # false_negatives = gm_annotations - true_positives
    false_negatives = gm_annotations
    for tp in true_positives:
        false_negatives = false_negatives.exclude(start=tp.start, text=tp.text)

    score = determine_f(len(true_positives), false_positives.count(), false_negatives.count())
    return (score, true_positives, false_positives, false_negatives)


def create_from_pubmed_id(pubmed_id=None):
    pubmed_id = str(pubmed_id)

    # Check if the account already exists
    try:
        return Document.objects.get(document_id=pubmed_id)

    except ObjectDoesNotExist:
        Entrez.email = settings.ENTREZ_EMAIL
        h = Entrez.efetch(db='pubmed', id=[pubmed_id], rettype='medline', retmode='text')
        records = Medline.parse(h)

        for record in records:
            # http://www.nlm.nih.gov/bsd/mms/medlineelements.html

            if record.get('TI') and record.get('PMID') and record.get('CRDT'):
                doc = Document()

                doc.document_id = record.get('PMID')
                doc.title = record.get('TI')
                doc.created = datetime.datetime.strptime(record.get('CRDT')[0], '%Y/%m/%d %H:%M')
                doc.source = 'pubmed'
                doc.save()

                sec = Section(kind='o')
                sec.document = doc
                sec.save()

                sec = Section(kind='t')
                sec.text = record.get('TI')
                sec.document = doc
                sec.save()

                if record.get('AB'):
                    sec = Section(kind='a')
                    sec.text = record.get('AB')
                    sec.document = doc
                    sec.save()

                return doc
            break


