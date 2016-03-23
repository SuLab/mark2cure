from django.db.models.signals import post_save
from django.dispatch import receiver

from brabeion import badges

@receiver(post_save, sender='score.Point')
def point_post_save(sender, instance, created, **kwargs):
    point = instance

    if created and not kwargs.get('raw', False):
        badges.possibly_award_badge("points_awarded", user=point.user)

