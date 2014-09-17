from django.conf import settings
from django.dispatch import receiver
from django.db.models.signals import post_save

from mark2cure.account.models import UserProfile

from django.core.mail import send_mail
import sys
import logging
logger = logging.getLogger(__name__)


@receiver(post_save, sender=UserProfile)
def message_post_save(sender, instance, **kwargs):
    userprofile = instance

    if 'test' not in sys.argv and settings.DEBUG is False:
        try:
            send_mail(
                '[Mark2Cure] New User Signup',
                'Email: {email_address}'.format(
                    email_address=userprofile.user.email),
                settings.SERVER_EMAIL,
                [email[1] for email in settings.MANAGERS])
        except:
            pass
