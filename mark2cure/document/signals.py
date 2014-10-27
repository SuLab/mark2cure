from django.conf import settings
from django.db.models import signals
from django.core.mail import send_mail

from mark2cure.document.models import Comment

import logging
logger = logging.getLogger(__name__)


def comment_save_handler(sender, instance, **kwargs):
    comment = instance

    send_mail('[Mark2Cure #{0}] Document comment'.format(settings.EXPERIMENT),
              '{0} commented: {1} on document id {2}'.format(comment.user.pk, comment.message, comment.document.pk),
              settings.SERVER_EMAIL,
              [email[1] for email in settings.MANAGERS])

signals.post_save.connect(comment_save_handler, sender=Comment)
