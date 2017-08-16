from django.contrib.auth.models import User
from django.contrib.contenttypes.models import ContentType

from ...analysis.tasks import generate_reports
from ...analysis.models import Report
from ...document.models import View, Annotation
from .models import EntityRecognitionAnnotation

import random


def select_best_opponent(task, document, player):
    '''
        Select the best opponate for a certain task meaning a
        certain document scoped to a quest for a certain task type.

        1) First weight by GM, if one is available always prefer it over
            other users

        2) Select users with non-empty response for this View (Document scoped to Quest)
            Explanation: This ensures we only look at users who have submitted the document
                         so that a comparison can be shown

            2.1) Select which available comparsion has the best F Score for this Group
                  Explanation: A user performance varies within Groups as they require
                               different types of skills

            Sort by internal F score average, select top 3 (over threshold), pick 1 at random
                * Sort by top weighted F-Scores

        3) If no GM or no non-empty responses return None

    Args:

    Returns:
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
    for quest_relationship in others_quest_relationships.exclude(user=gm_user).exclude(user__pk__in=[107, ]):
        if quest_relationship.views.filter(section__document=document, completed=True).exists():
            # (TODO) Don't add option of them unless they've submitted 1+ Annotations
            view_ids = quest_relationship.views.filter(section__document=document, completed=True).values_list('pk', flat=True)
            if Annotation.objects.filter(kind='e', view__pk__in=view_ids).exists():
                previous_users.append(quest_relationship.user)

    if others_quest_relationships.exists() and len(previous_users):
        # No expert around so select a previous user at random
        previous_users_pks = [str(u.pk) for u in previous_users]

        report = task.group.report_set.filter(report_type=Report.AVERAGE).order_by('-created').first()
        if report:

            # (TODO) if str load up manually
            df = report.dataframe
            df = df[df['user_id'].isin(previous_users_pks)]

            if df.shape[0] == 0:
                previous_users_pks = [int(u.pk) for u in previous_users]
                df = report.dataframe
                df = df[df['user_id'].isin(previous_users_pks)]

            row_length = df.shape[0]
            if row_length:
                # Top 1/2 of the users (sorted by F)
                df = df.iloc[:int(row_length / 2)]

                # Select 1 at random
                top_half_user_pks = list(df.user_id)
                random.shuffle(top_half_user_pks)
                selected_user_pk = top_half_user_pks[0]

                for u in previous_users:
                    if str(u.pk) == str(selected_user_pk):
                        return u
        else:
            generate_reports.apply_async(
                args=[task.group.pk],
                queue='mark2cure_tasks')

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
    content_type_id = str(ContentType.objects.get_for_model(EntityRecognitionAnnotation.objects.all().first()).id)
    gm_annotations = EntityRecognitionAnnotation.objects.annotations_for_view_pks([v.pk for v in gm_views if type(v) is View], content_type_id)
    user_annotations = EntityRecognitionAnnotation.objects.annotations_for_view_pks([v.pk for v in user_views if type(v) is View], content_type_id)

    true_positives = [gm_ann for gm_ann in gm_annotations if match_exact(gm_ann, user_annotations)]

    # print 'True Pos', len(true_positives)
    # print 'User Anns', len(user_annotations)

    # Annotations the user submitted that were wrong (the User set without their True Positives)
    # false_positives = user_annotations - true_positives
    false_positives = user_annotations
    for tp in true_positives:
        false_positives = list(filter(lambda er_ann: er_ann.start != tp.start and er_ann.text != tp.text, false_positives))
    # print 'False Positives', len(false_positives)

    # # Annotations the user missed (the GM set without their True Positives)
    # false_negatives = gm_annotations - true_positives
    # (TODO) FN appears to be 1/2 what it should in most cases --Max 3/23/16
    false_negatives = gm_annotations
    for tp in true_positives:
        false_negatives = list(filter(lambda er_ann: er_ann.start != tp.start and er_ann.text != tp.text, false_negatives))
    # print 'False Negatives', len(false_negatives)

    score = determine_f(len(true_positives), len(false_positives), len(false_negatives))
    return (score, true_positives, false_positives, false_negatives)
