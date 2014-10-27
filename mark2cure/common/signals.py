from django.dispatch import receiver
from django.db.models.signals import post_save
from django.conf import settings

from mark2cure.common.models import UserQuestRelationship

from brabeion import badges
from brabeion.models import BadgeAward

import sys
import os


@receiver(post_save, sender=UserQuestRelationship)
def user_quest_relationship_post_save(sender, instance, **kwargs):
    user_quest_rel = instance

    print "TRIGGERED", user_quest_rel
    #user_quest_rel.completed = True
    #user_quest_rel.save()

    #if user_quest_rel.completed:
    #    user_quest_rel.user.profile.rating.add(score=user_quest_rel.task.points, user=None, ip_address=os.urandom(7).encode('hex'))
    #    badges.possibly_award_badge("points_awarded", user=user_quest_rel.user)


