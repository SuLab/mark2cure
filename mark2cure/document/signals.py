from django.conf import settings
from django.db.models import signals
from django.core.mail import send_mail

from mark2cure.document.models import Document, Activity, Comment

import requests, importlib, logging
logger = logging.getLogger(__name__)


def comment_save_handler(sender, instance, **kwargs):
    comment = instance

    send_mail('[Mark2Cure #{0}] Document comment',
              '{1} commented: {2} on document id {3}'.format(settings.EXPERIMENT, comment.user.pk, comment.message, comment.document.pk),
              settings.SERVER_EMAIL,
              [email[1] for email in settings.MANAGERS])


def activity_save_handler(sender, instance, **kwargs):
    activity = instance
    user = activity.user
    user_profile = user.userprofile
    '''
      Email notify monitors
    '''
    if activity.f_score == 1.0 or activity.f_score == 0.0:
        send_mail('[Mark2Cure #{0}] HIT completion',
                      '{1} scored {2} on document id {3}'.format(settings.EXPERIMENT, user.pk, activity.f_score, activity.document.pk),
                       settings.SERVER_EMAIL,
                       [email[1] for email in settings.MANAGERS])


    '''
      Softblock on poor performance
    '''
    if user_profile.mturk and activity.submission_type == "gm":
        if activity.f_score <= 0.5:

            if user_profile.mturk:
                latest_results = Activity.objects.filter(user=user, experiment=settings.EXPERIMENT)[:3]
            else:
                latest_results = Activity.objects.filter(user=user)[:3]

            if len(latest_results) == 3:
                if poor_subs_count[1].f_score < .5 and poor_subs_count[2].f_score < .5:
                    user_profile.softblock = True
                    user_profile.save()
                    send_mail('[Mark2Cure #{0}] softblock',
                                '{1} was blocked due to document id {3}'.format(settings.EXPERIMENT, user.pk, activity.document.pk),
                                settings.SERVER_EMAIL,
                                [email[1] for email in settings.MANAGERS])

signals.post_save.connect(activity_save_handler, sender=Activity)
signals.post_save.connect(comment_save_handler, sender=Comment)


