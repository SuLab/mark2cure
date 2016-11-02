from django.db.models.signals import post_save
from django.dispatch import receiver
from .tasks import submit_pubtator


@receiver(post_save, sender='document.Pubtator')
def pubtator_post_save(sender, instance, created, **kwargs):
    '''
        When a new Pubtator model is created,
        start fetching it's content
    '''
    pubtator = instance

    if created is True and pubtator.content is None and not kwargs.get('raw', False):
        submit_pubtator.apply_async(
            args=[instance.pk],
            queue='mark2cure_tasks')

