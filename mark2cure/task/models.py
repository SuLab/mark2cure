from django.contrib.auth.models import User
from django.db import models


class Level(models.Model):
    user = models.ForeignKey(User)

    ENTITY_RECOGNITION = 'e'
    RELATION = 'r'
    TASK_TYPE_CHOICES = (
        (ENTITY_RECOGNITION, 'Entity Recognition'),
        (RELATION, 'Relation'),
    )
    task_type = models.CharField(max_length=1, choices=TASK_TYPE_CHOICES)
    level = models.IntegerField()
    created = models.DateTimeField(auto_now_add=True)

    class Meta:
        get_latest_by = 'created'
        ordering = ('-level', )

    def get_name(self):
        if self.task_type == 'e':
            levels = ["Basic", "Disease Marking", "Disease Advanced", "Disease Matching", "Intermediate", "Proficient", "Advanced", "Expert"]
        else:
            levels = ["beginner", "Medium", "Expert"]

        return levels[self.level]


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
    documents = models.ManyToManyField('document.Document', through='DocumentQuestRelationship', blank=True)
    users = models.ManyToManyField(User, through='UserQuestRelationship', blank=True)
    points = models.IntegerField(max_length=6, blank=True, default=0)
    experiment = models.IntegerField(blank=True, null=True)

    requires_qualification = models.IntegerField(max_length=6, blank=True, null=True)
    provides_qualification = models.IntegerField(max_length=6, blank=True, null=True)
    meta_url = models.CharField(max_length=200, null=True, blank=True)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    # Tasks are not shared between groups so no need for m2m
    group = models.ForeignKey('common.Group', blank=True, null=True)

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
                from ..document.models import View
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

    views = models.ManyToManyField('document.View')

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
        return u'/task/entity-recognition/quest/{quest_pk}/ {username}'.format(
            quest_pk=self.task.pk,
            username=self.user.username)


class DocumentQuestRelationship(models.Model):
    task = models.ForeignKey(Task)
    document = models.ForeignKey('document.Document')

    points = models.IntegerField(max_length=7, blank=True, default=0)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        return u'Document Quest Relationship'

