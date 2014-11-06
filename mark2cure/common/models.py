from django.db import models
from django.conf import settings
from django.core.validators import MaxValueValidator, MinValueValidator

from mark2cure.document.models import Document, View
from django.contrib.auth.models import User

from brabeion import badges
from brabeion.base import Badge, BadgeAwarded


class SkillBadge(Badge):
    slug = "skill"
    levels = [
        "Basic",
        "Disease Marking"
        "Disease Advanced"
        "Intermediate",
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
        if level and level == current_highest+2:
            return BadgeAwarded(level=level)


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


class Task(models.Model):
    name = models.CharField(max_length=200)

    TRAINING = 't'
    QUEST = 'q'
    KIND_CHOICES = (
        (TRAINING, 'Training'),
        (QUEST, 'Quest'),
    )
    kind = models.CharField(max_length=1, choices=KIND_CHOICES, default=QUEST)

    completions = models.IntegerField(default=10)
    documents = models.ManyToManyField(Document, through='DocumentQuestRelationship', blank=True)
    users = models.ManyToManyField(User, through='UserQuestRelationship', blank=True)
    points = models.IntegerField(max_length=6, blank=True, default=50)

    requires_qualification = models.IntegerField(max_length=6, blank=True, null=True)
    provides_qualification = models.IntegerField(max_length=6, blank=True, null=True)
    meta_url = models.CharField(max_length=200, null=True, blank=True)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    def total_points(self):
        print self.points + sum(DocumentQuestRelationship.objects.filter(task=self).values_list('points', flat=True))

    def create_views(self, document, user):
        user_quest_rel_views = self.userquestrelationship_set.get(user=user).views

        #if user_quest_rel_views.filter(section__document=document).count() < document.count_available_sections():
        if user_quest_rel_views.filter(section__document=document).count() < document.count_available_sections():

            for sec in document.available_sections():
                view = View.objects.create(section=sec, user=user)
                user_quest_rel_views.add(view)


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

    def __unicode__(self):
        return u'User Quest Relationship'

class DocumentQuestRelationship(models.Model):
    task = models.ForeignKey(Task)
    document = models.ForeignKey(Document)

    points = models.IntegerField(max_length=7, blank=True, default=5)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        return u'Document Quest Relationship'

