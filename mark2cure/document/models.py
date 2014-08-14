from django.db import models
from django.conf import settings
from django.core.validators import MaxValueValidator, MinValueValidator
from django.contrib.auth.models import User

from mark2cure.document.managers import DocumentManager

from decimal import Decimal as D
from nltk.tokenize import WhitespaceTokenizer

import random, itertools


class Document(models.Model):
    document_id = models.IntegerField(blank=True)
    title = models.TextField(blank=False)
    authors = models.TextField(blank=False)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)
    source = models.CharField(max_length=200, blank=True)

    objects = DocumentManager()

    def __unicode__(self):
        return self.title

    def submitted(self):
        return Activity.objects.filter(document=self).count()

    def available_sections(self):
        return self.section_set.exclude(kind='o').all()


    def count_available_sections(self):
        return self.section_set.exclude(kind='o').count()

    def has_golden(self):
        # (TODO) Change so that documents can be golden or normal cases, don't tie to the exisitance of a user's annotations
        return Annotation.objects.filter(view__user__username='goldenmaster', view__section__document = self).exists()


    def is_complete(self, user, user_profile, sections, task_type='cr'):
        # Stick w/ Views b/c the Activity results haven't been logged yet
        query = View.objects.filter(user=user, completed=True, task_type=task_type, section__document=self, experiment=settings.EXPERIMENT if user_profile.mturk else None)
        return True if query.count() >= len(sections) else False


    def create_views(self, user, task_type, completed=False):
        '''
          We never want to 'recycle' a View, always make a new one, even if there
          are currently 'open' / uncompleted ones from the same doc
        '''
        for sec in self.available_sections():
            view = View(task_type=task_type, section=sec, user=user, experiment=settings.EXPERIMENT if user.userprofile.mturk else None)
            view.completed = completed
            view.save()


    def update_views(self, user, task_type, completed=False):
        for sec in self.available_sections():
            view = View.objects.filter(user=user, task_type=task_type, section=sec).latest()
            view.completed = completed
            view.save()

    def latest_views(self, user, task_type='cr', completed=True):
        return View.objects.filter(user=user, task_type=task_type, completed=completed, experiment=settings.EXPERIMENT if user.userprofile.mturk else None, section__document=self)[:2]

    def latest_annotations(self, user=None):
        if user is None:
            user = User.objects.get(username='goldenmaster')
            return Annotation.objects.filter(view__section__document=self, view__task_type='cr', kind='e', view__user=user).order_by('start')
        else:
            user_views = self.latest_views(user)
            return Annotation.objects.filter(view__pk__in=[view.pk for view in user_views]).order_by('start')

    class Meta:
        ordering = ('-created',)
        get_latest_by = 'updated'


class Section(models.Model):
    SECTION_KIND_CHOICE = (
      ('o', 'Overview'),
      ('t', 'Title'),
      ('a', 'Abstract'),
      ('p', 'Paragraph'),
      ('f', 'Figure'),
    )
    kind = models.CharField(max_length=1, choices=SECTION_KIND_CHOICE)
    text = models.TextField(blank=True)
    source = models.ImageField(blank=True, upload_to='media/images/', default='images/figure.jpg')
    cache = models.TextField(blank=True)
    concepts = models.ManyToManyField('Concept', blank=True, null=True)
    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    document = models.ForeignKey(Document)


    def latest_view(self, user, task_type='cr', completed=True):
        return View.objects.filter(user=user, task_type=task_type, completed=completed, experiment=settings.EXPERIMENT if user.userprofile.mturk else None, section=self).first()


    def latest_annotations(self, user = None):
        if user is None:
            user = User.objects.get(username = 'goldenmaster')
            return Annotation.objects.filter(view__section=self, view__task_type='cr', kind='e', view__user=user).order_by('start')
        else:
            # (TODO) Scope to section only
            user_view = self.latest_view(user)
            return Annotation.objects.filter(view__pk = user_view.pk).order_by('start')


    def resultwords(self, user):
        # Gather words and positions from the text
        words_index = WhitespaceTokenizer().span_tokenize(self.text)
        words_text = WhitespaceTokenizer().tokenize(self.text)
        words = zip(words_index, words_text)

        # Add counters for concensus count and personal annotation
        # ((start, stop), 'Word string itself', Intensity, GM Ann ID, User Ann ID, Did user annotate)
        words = [w + (0, None, None, False,) for w in words]

        # Gather other annotations from GM and users for this section
        gm_anns = self.latest_annotations().values_list('pk', 'start', 'text')
        # Build the running counter of times a word was annotated
        for gm_pk, start, text in gm_anns:
          length = len(text)

          for idx, word in enumerate(words):
            word_start = word[0][0]
            counter = word[2]
            if word_start >= start and word_start <= start + length:
              counter += 1
              words[idx] = (word[0], word[1], counter, gm_pk, word[3], word[4])


        user_anns = self.latest_annotations(user).values_list('pk', 'start', 'text')
        # Build the running counter of times a word was annotated
        for user_pk, start, text in user_anns:
          length = len(text)
          for idx, word in enumerate(words):
            word_start = word[0][0]
            if word_start >= start and word_start <= start + length:
              words[idx] = (word[0], word[1], word[2], word[3], user_pk, True)

        return words

    def update_view(self, user, task_type, completed=False):
        view = View.objects.filter(user=user, task_type=task_type, section=self).latest()
        view.completed = completed
        view.save()
        return view

    def __unicode__(self):
        if self.kind == 'o':
          return u'[Overview] {0}'.format(self.document.title)
        else:
          return self.text

    class Meta:
        get_latest_by = 'updated'


TASK_TYPE_CHOICE = (
  ('cr', 'Concept Recognition'),
  ('cn', 'Concept Normalization'),
  ('rv', 'Relationship Verification'),
  ('ri', 'Relationship Identification'),
  ('rc', 'Relationship Correction'),
)

class View(models.Model):
    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    task_type = models.CharField(max_length=3, choices=TASK_TYPE_CHOICE, blank=True, default='cr')
    completed = models.BooleanField(default=False, blank=True)
    experiment = models.IntegerField(blank=True, null=True)

    section = models.ForeignKey(Section)
    user = models.ForeignKey(User)

    def __unicode__(self):
        return u'Doc:'+ str(self.section.document.pk) +', Sec:'+ str(self.section.pk) +' by '+ self.user.username

    class Meta:
        get_latest_by = 'updated'
        ordering = ('-updated',)


class Activity(models.Model):
    '''
      The way to store out progress of completions and scores of worker permformance
    '''
    user = models.ForeignKey(User)
    document = models.ForeignKey(Document)
    task_type = models.CharField(max_length=3, choices=TASK_TYPE_CHOICE, blank=True, default='cr')
    experiment = models.IntegerField(blank=True, null=True)

    SUBMISSION_TYPE = (
      ('gm', 'Golden Master'),
      ('cc', 'Community Consensus'),
      ('na', 'Never Annotated'),
    )
    submission_type = models.CharField(max_length=3, choices=SUBMISSION_TYPE, blank=True, default='gm')

    precsion = models.DecimalField(max_digits=11, decimal_places=5, validators=[MaxValueValidator(1), MinValueValidator(0)], null=True, blank=True)
    recall = models.DecimalField(max_digits=11, decimal_places=5, validators=[MaxValueValidator(1), MinValueValidator(0)], null=True, blank=True)
    f_score = models.DecimalField(max_digits=11, decimal_places=5, validators=[MaxValueValidator(1), MinValueValidator(0)], null=True, blank=True)

    created = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ('-created',)
        get_latest_by = 'created'

    def __unicode__(self):
        return u'Activity {0}'.format(self.user)


class Refute(models.Model):
    '''
      I ended up putting this into a separate model after careful thought:
       - There may be multiple types of refutes in the future
       - Simple query for # Refutes per section
       - Separate update/create timestamps, makes a difference if we add
       a `resolved` boolean in the future
    '''
    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)
    message = models.TextField(blank=True)
    view = models.ForeignKey(View)

    def __unicode__(self):
        return '{0} {1}'.format(self.message, self.view)

    class Meta:
        get_latest_by = 'updated'


class Comment(models.Model):
    '''
      Refutes without dedicated sections
    '''
    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    message = models.TextField(blank=True)

    task_type = models.CharField(max_length=3, choices=TASK_TYPE_CHOICE, blank=True, default='cr')
    document = models.ForeignKey(Document)
    user = models.ForeignKey(User)

    def __unicode__(self):
        return u'{0} {1} {2}'.format(self.message, self.document, self.user)

    class Meta:
        get_latest_by = 'updated'


class Annotation(models.Model):
    ANNOTATION_KIND_CHOICE = (
      ('e', 'Entities'),
      ('a', 'Attributes'),
      ('r', 'Relations'),
      ('t', 'Triggers'),
      ('e', 'Events'),
    )
    kind = models.CharField(max_length=1, choices=ANNOTATION_KIND_CHOICE, blank=False)

    # Disease, Gene, Protein, et cetera...
    type = models.CharField(max_length=40, blank=True, null=True, default='disease')

    text = models.TextField(blank=True, null=True)
    start = models.IntegerField(blank=True, null=True)

    created = models.DateTimeField(auto_now_add=True)

    view = models.ForeignKey(View)

    def is_exact_match(self, comparing_annotation):
        return True if self.start == comparing_annotation.start and len(self.text) == len(comparing_annotation.text) else False

    def __unicode__(self):
        if self.kind == 'r':
            return 'Relationship Ann'
        if self.kind == 'e':
            return '{0} ({1}) [{2}]'.format(self.text, self.start, self.pk)
        else:
            return self.type

    class Meta:
        get_latest_by = 'updated'
        # (TODO) This is not supported by MySQL but would help prevent dups in this table
        # unique_together = ['kind', 'type', 'text', 'start', 'view']


class Concept(models.Model):
    concept_id = models.TextField(blank=False)
    preferred_name = models.TextField(blank=True)
    definition = models.TextField(blank=True)
    semantic_type = models.TextField(blank=True)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        return '{0} ({1})'.format(self.preferred_name, self.concept_id)


class RelationshipType(models.Model):
    full_name = models.CharField(max_length=80)
    type = models.CharField(max_length=80)
    definition = models.TextField(blank=True)

    parent = models.ForeignKey('self', blank=True, null=True, related_name='children')

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        names = []
        current = self.parent
        while current:
            names.append(current.full_name)
            current = current.parent

        names = list(reversed(names))
        path = ' :: '.join(names)

        return '{0}:: {1}'.format(path, self.full_name)


class ConceptRelationship(models.Model):
    concept = models.ForeignKey(Concept, related_name='actor')
    relationship = models.ForeignKey('RelationshipType', blank=True, null=True)
    target = models.ForeignKey(Concept, blank=True, null=True, related_name='target')

    annotation = models.ForeignKey(Annotation, blank=True, null=True)

    validate = models.ForeignKey('self', blank=True, null=True)
    confidence = models.DecimalField(max_digits=11, decimal_places=5, validators=[MaxValueValidator(1), MinValueValidator(0)], default=0)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        return '{0} >> {1} >> {2}'.format(self.concept.preferred_name,
                                          self.relationship.full_name,
                                          self.target)


