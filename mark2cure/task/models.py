from django.contrib.auth.models import User
from django.db import models


ENTITY_RECOGNITION = 'ner'
RELATION = 're'
TASK_TYPE_CHOICES = (
    (ENTITY_RECOGNITION, 'Name Entity Recognition'),
    (RELATION, 'Relation Extraction'),
)


class Requirement(models.Model):
    """Required levels (identifiers) for a task to be enabled
    """
    # (TODO) why can't I call this id???
    hash = models.CharField(max_length=8)

    task_type = models.CharField(max_length=3, choices=TASK_TYPE_CHOICES)
    name = models.CharField(max_length=200)
    order = models.IntegerField()

    # For training levels no longer used
    active = models.BooleanField(default=False)

    class Meta:
        get_latest_by = 'order'
        ordering = ('-order', )


class Level(models.Model):
    """The Task Type specific Level a user is trained at
    """
    user = models.ForeignKey(User)
    requirement = models.ForeignKey('task.Requirement', blank=True, null=True, related_name="completes")
    created = models.DateTimeField(auto_now_add=True)

    class Meta:
        get_latest_by = 'created'
        ordering = ('-created', )

    def get_name(self):
        if self.requirement.task_type == 'ner':
            levels = ["Basic", "Disease Marking", "Disease Advanced", "Disease Matching", "Intermediate", "Proficient", "Advanced", "Expert", "Expert", "Expert"]
        else:
            levels = ["Beginner", "Medium", "Expert", "Expert", "Expert"]

        return levels[self.requirement.order]


class Task(models.Model):
    """This is an ER Quest, tracks whose completed it and what documents are contained within
        the Quest

        * Originally called Task when training was going to be dynamic, a TODO is remove all
            'Training' references from this model
    """
    name = models.CharField(max_length=200)

    # If no completions defined, allow infinity K value
    completions = models.IntegerField(default=10, blank=True, null=True)
    documents = models.ManyToManyField('document.Document', through='DocumentQuestRelationship', blank=True)
    users = models.ManyToManyField(User, through='UserQuestRelationship', blank=True)
    points = models.IntegerField(blank=True, default=0)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    # Tasks are not shared between groups so no need for m2m
    group = models.ForeignKey('common.Group', blank=True, null=True)

    def create_views(self, document, user):
        user_quest_rel = self.userquestrelationship_set.filter(user=user, completed=False).first()
        user_quest_rel_views = user_quest_rel.views
        if user_quest_rel_views.filter(section__document=document).count() < document.count_available_sections():
            for sec in document.available_sections():
                from ..document.models import View
                view = View.objects.create(section=sec, user=user)
                user_quest_rel_views.add(view)

    def __unicode__(self):
        return self.name


class UserQuestRelationship(models.Model):
    """Describes a User's Status on a specific ER Quest

        We use this to track if a User has completed a Quest
        * Technically this can be inferred by if they'ved completed all of the
          documents within the quest, but was originally made to separate if a
          PMID appears in multiple Quests
    """
    task = models.ForeignKey(Task)
    user = models.ForeignKey(User)

    views = models.ManyToManyField('document.View', help_text='The 10 (5 Quests * 2 Sections) Viewable blocks of text that are part of the Quest')

    completed = models.BooleanField(default=False, blank=True)
    score = models.IntegerField(blank=True, default=5)

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

    points = models.IntegerField(blank=True, default=0)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        return u'Document Quest Relationship'

