from django.core.management.base import BaseCommand
from django.contrib.auth.models import User
from mark2cure.common.models import Task, UserQuestRelationship
from brabeion import badges
from mark2cure.score.models import Point

import os


class Command(BaseCommand):
    help = 'Assign the Username all the required training tests'

    def handle(self, *args, **options):
        user = User.objects.get(username=args[0])
        print user


        badges.possibly_award_badge("skill_awarded", user=user, level=1)

        task = Task.objects.first()
        # Required importing
        from django.contrib.contenttypes.models import ContentType
        from django.utils import timezone
        content_type = ContentType.objects.get_for_model(task)

        Point.objects.create(user=user, amount=task.points, content_type=content_type, object_id=task.id, created=timezone.now())
        UserQuestRelationship.objects.create(task=task, user=user, completed=True)
        badges.possibly_award_badge("skill_awarded", user=user, level=2)

        task = Task.objects.get(pk=2)
        UserQuestRelationship.objects.create(task=task, user=user, completed=True)
        Point.objects.create(user=user, amount=task.points, content_type=content_type, object_id=task.id, created=timezone.now())
        badges.possibly_award_badge("skill_awarded", user=user, level=3)

        task = Task.objects.get(pk=3)
        UserQuestRelationship.objects.create(task=task, user=user, completed=True)
        Point.objects.create(user=user, amount=task.points, content_type=content_type, object_id=task.id, created=timezone.now())
        badges.possibly_award_badge("skill_awarded", user=user, level=4)

