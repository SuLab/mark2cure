from django.core.management.base import BaseCommand
from django.contrib.auth.models import User
from mark2cure.common.models import Task, UserQuestRelationship
from mark2cure.score.models import Point

from mark2cure.task.models import Level


class Command(BaseCommand):
    help = 'Assign the Username all the required training tests'

    def handle(self, *args, **options):
        user = User.objects.get(username=args[0])

        # Required importing
        from django.contrib.contenttypes.models import ContentType
        from django.utils import timezone

        Level.objects.create(user=user, task_type='ner', level=1, created=timezone.now())

        task = Task.objects.first()
        content_type = ContentType.objects.get_for_model(task)

        Point.objects.create(user=user, amount=task.points, content_type=content_type, object_id=task.id, created=timezone.now())
        UserQuestRelationship.objects.create(task=task, user=user, completed=True)
        Level.objects.create(user=user, task_type='ner', level=2, created=timezone.now())

        task = Task.objects.get(pk=2)
        UserQuestRelationship.objects.create(task=task, user=user, completed=True)
        Point.objects.create(user=user, amount=task.points, content_type=content_type, object_id=task.id, created=timezone.now())
        Level.objects.create(user=user, task_type='ner', level=3, created=timezone.now())

        task = Task.objects.get(pk=3)
        UserQuestRelationship.objects.create(task=task, user=user, completed=True)
        Point.objects.create(user=user, amount=task.points, content_type=content_type, object_id=task.id, created=timezone.now())
        Level.objects.create(user=user, task_type='ner', level=4, created=timezone.now())

