from django.core.management.base import BaseCommand
from django.contrib.auth.models import User
from mark2cure.common.models import Task, UserQuestRelationship
from brabeion import badges

import os


class Command(BaseCommand):
    help = 'Assign the Username all the required training tests'

    def handle(self, *args, **options):
        user = User.objects.get(username=args[0])
        print user

        badges.possibly_award_badge("skill_awarded", user=user, level=1)

        task = Task.objects.first()
        user.profile.rating.add(score=task.points, user=None, ip_address=os.urandom(7).encode('hex'))
        UserQuestRelationship.objects.create(task=task, user=user, completed=True)
        badges.possibly_award_badge("skill_awarded", user=user, level=2)

        task = Task.objects.get(pk=2)
        UserQuestRelationship.objects.create(task=task, user=user, completed=True)
        user.profile.rating.add(score=task.points, user=None, ip_address=os.urandom(7).encode('hex'))
        badges.possibly_award_badge("points_awarded", user=user)
        badges.possibly_award_badge("skill_awarded", user=user, level=3)

        task = Task.objects.get(pk=3)
        UserQuestRelationship.objects.create(task=task, user=user, completed=True)
        user.profile.rating.add(score=task.points, user=None, ip_address=os.urandom(7).encode('hex'))
        badges.possibly_award_badge("points_awarded", user=user)
        badges.possibly_award_badge("skill_awarded", user=user, level=4)

