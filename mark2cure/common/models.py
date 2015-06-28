from django.db import models
from django.db.models import Count

from ..document.models import Document, View
from django.contrib.auth.models import User

from brabeion import badges
from brabeion.base import Badge, BadgeAwarded

from decimal import Decimal


class SkillBadge(Badge):
    slug = "skill"
    levels = [
        "Basic",
        "Disease Marking",  # T1 complete
        "Disease Advanced",  # T2 complete
        "Disease Matching",  # T3 complete
        "Intermediate",  # 1st GM Quest Complete
        "Proficient",
        "Advanced",
        "Expert",
    ]
    events = [
        "skill_awarded",
    ]
    multiple = False

    def award(self, **state):
        user = state["user"]
        level = state.pop("level", None)
        current_highest = user.profile.highest_level(self.slug).level

        if (level and level == current_highest + 1) or state.get('force', None):
            return BadgeAwarded(level=level + 1)

badges.register(SkillBadge)


class PointsBadge(Badge):
    slug = "points"
    levels = [
        "Novice Magic Marker",
        "Magic Marker",
        "Expert Magic Marker",
        "Novice Marker Bee",
        "Marker Bee",
        "Expert Marker Bee",
        "Novice Marksman",
        "Marksman",
        "Expert Marksman",
        "Novice Mark up",
        "Mark up",
        "Expert Mark up"
        "Novice Benchmarker",
        "Benchmarker",
        "Expert Benchmarker",
        "Master Benchmarker",
        "Master Marker",
        "Bronze Master Marker",
        "Silver Master Marker",
        "Gold Master Marker",
    ]
    events = [
        "points_awarded",
    ]
    multiple = False

    def award(self, **state):
        user = state["user"]
        points = user.profile.rating_score
        if points > 500:
            return BadgeAwarded(level=2)
        if points > 50:
            return BadgeAwarded(level=1)


badges.register(PointsBadge)


class Group(models.Model):
    name = models.CharField(max_length=200)
    stub = models.CharField(max_length=200)
    description = models.TextField(blank=True, null=True)
    order = models.DecimalField(default=0, max_digits=3, decimal_places=3)

    enabled = models.BooleanField(default=False)

    def get_documents(self):
        # (TODO?) Return for __in of task_ids
        return Document.objects.filter(task__group=self)

    def percentage_complete(self):
        task_queryset = self.task_set.extra(select = {
        "completed" : """
            SELECT COUNT(*) AS completed
            FROM common_userquestrelationship
            WHERE (common_userquestrelationship.completed = 1
                AND common_userquestrelationship.task_id = common_task.id)"""
            })
        completed = sum(task_queryset.values_list('completed', flat=True))
        required = sum(task_queryset.values_list('completions', flat=True))
        if required:
            return (Decimal(completed) / Decimal(required))*100
        else:
            return 0

    def __unicode__(self):
        return self.name


class Task(models.Model):
    name = models.CharField(max_length=200)

    TRAINING = 't'
    QUEST = 'q'
    KIND_CHOICES = (
        (TRAINING, 'Training'),
        (QUEST, 'Quest'),
    )
    kind = models.CharField(max_length=1, choices=KIND_CHOICES, default=QUEST)

    # If no completions defined, allow infinity K value
    completions = models.IntegerField(default=10, blank=True, null=True)
    documents = models.ManyToManyField(Document, through='DocumentQuestRelationship', blank=True)
    users = models.ManyToManyField(User, through='UserQuestRelationship', blank=True)
    points = models.IntegerField(max_length=6, blank=True, default=0)
    experiment = models.IntegerField(blank=True, null=True)

    requires_qualification = models.IntegerField(max_length=6, blank=True, null=True)
    provides_qualification = models.IntegerField(max_length=6, blank=True, null=True)
    meta_url = models.CharField(max_length=200, null=True, blank=True)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    # Tasks are not shared between groups so no need for m2m
    group = models.ForeignKey(Group, blank=True, null=True)

    def total_points(self):
        print self.points + sum(DocumentQuestRelationship.objects.filter(task=self).values_list('points', flat=True))

    def remaining_documents_count(self, completed_ids=[]):
        # Document pks for the Quest the user has left
        if self.documents:
            # Document pks for the Quest the user has left
            return self.documents.exclude(pk__in=completed_ids).count()
        else:
            return 0

    def remaining_documents(self, completed_ids=[]):
        # Document pks for the Quest the user has left
        if self.documents:
            # Document pks for the Quest the user has left
            return list(self.documents.exclude(pk__in=completed_ids).all())
        else:
            return []


    def user_relationship(self, user, completed=False):
        return UserQuestRelationship.objects.filter(task=self, user=user, completed=completed).first()


    def create_views(self, document, user):
        user_quest_rel = self.userquestrelationship_set.filter(user=user, completed=False).first()
        user_quest_rel_views = user_quest_rel.views
        if user_quest_rel_views.filter(section__document=document).count() < document.count_available_sections():

            for sec in document.available_sections():
                view = View.objects.create(section=sec, user=user)
                user_quest_rel_views.add(view)

    def complete_views(self, document, user):
        user_quest_rel = self.userquestrelationship_set.filter(user=user, completed=False).first()
        user_quest_rel_views = user_quest_rel.views

        for view in user_quest_rel_views.filter(section__document=document).all():
            # (TODO) Validate (require 1+ ann for example?) before allowing completion
            view.completed = True
            view.save()

    def clear_documents(self):
        # Remove any previous DocumentQuestRelationship which may have been in place
        for doc in self.documents.all():
            dqr = DocumentQuestRelationship.objects.get(document=doc, task=self)
            dqr.delete()

    def __unicode__(self):
        return self.name


class UserQuestRelationship(models.Model):
    task = models.ForeignKey(Task)
    user = models.ForeignKey(User)

    views = models.ManyToManyField(View)

    completed = models.BooleanField(default=False, blank=True)
    score = models.IntegerField(max_length=7, blank=True, default=5)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    def completed_views(self):
        return self.views.filter(completed=True)

    def completed_document_ids(self):
        # Collect the completed document pks the user has done for this quest
        return list(set(
            self.completed_views().values_list('section__document', flat=True)
        ))

    class Meta:
        get_latest_by = 'updated'


    def __unicode__(self):
        return u'/quest/{quest_pk}/ {username}'.format(
                quest_pk=self.task.pk,
                username=self.user.username)


class DocumentQuestRelationship(models.Model):
    task = models.ForeignKey(Task)
    document = models.ForeignKey(Document)

    points = models.IntegerField(max_length=7, blank=True, default=0)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        return u'Document Quest Relationship'


class SupportMessage(models.Model):
    user = models.ForeignKey(User, blank=True, null=True)
    text = models.TextField()
    referral = models.CharField(max_length=100, blank=True, null=True)
    created = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        return u'{text} (via {user})'.format(text=self.text, user=self.user)
