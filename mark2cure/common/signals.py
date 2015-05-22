from django.dispatch import receiver
from django.db.models.signals import post_save

from django.core.mail import EmailMessage

from .models import SupportMessage


@receiver(post_save, sender=SupportMessage)
def provider_post_save(sender, instance, **kwargs):
    message = instance

    msg = EmailMessage(to=['contact@mark2cure.org'])
    msg.template_name = 'mark2cure-support'
    msg.global_merge_vars = {
        'NAME': 'Unauthenicated User',
        'TEXT': message.text,
        'REFERRAL': message.referral
    }

    if message.user:
        msg.cc=[message.user.email]
        msg.global_merge_vars['NAME'] = message.user.username

        msg.send()
