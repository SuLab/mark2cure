from django.conf import settings
from django.db.models import signals
from django.core.mail import send_mail

from mark2cure.document.models import Activity, Comment

import logging
logger = logging.getLogger(__name__)


def comment_save_handler(sender, instance, **kwargs):
    comment = instance

    send_mail('[Mark2Cure #{0}] Document comment'.format(settings.EXPERIMENT),
              '{0} commented: {1} on document id {2}'.format(comment.user.pk, comment.message, comment.document.pk),
              settings.SERVER_EMAIL,
              [email[1] for email in settings.MANAGERS])


def activity_save_handler(sender, instance, **kwargs):
    activity = instance
    user = activity.user
    user_profile = user.userprofile

    '''
      Ban on poor performance
    '''
    if user_profile.mturk and activity.submission_type == 'gm':
        if activity.f_score <= 0.5:
            latest_results = Activity.objects.filter(user=user, experiment=settings.EXPERIMENT if user.userprofile.mturk else None, task_type='cr', submission_type='gm')[:3]

            if len(latest_results) == 3:
                if latest_results[1].f_score < .5 and latest_results[2].f_score < .5:
                    user_profile.banned = True
                    user_profile.save()
                    '''
                    send_mail('[Mark2Cure #{0}] banned'.format(settings.EXPERIMENT),
                                '{0} was blocked due to document id {1}'.format(user.pk, activity.document.pk),
                                settings.SERVER_EMAIL,
                                [email[1] for email in settings.MANAGERS])
                    '''


signals.post_save.connect(activity_save_handler, sender=Activity)
signals.post_save.connect(comment_save_handler, sender=Comment)
