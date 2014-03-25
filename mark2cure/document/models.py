from django.db import models
from django.conf import settings
from django.core.validators import MaxValueValidator, MinValueValidator
from django.utils.encoding import smart_text
from django.shortcuts import get_object_or_404

from mark2cure.document.managers import DocumentManager, AnnotationManager
from django.contrib.auth.models import User

from ttp import ttp
from decimal import Decimal as D
from copy import copy
from nltk.tokenize import WhitespaceTokenizer

import requests, random, datetime, itertools

class Document(models.Model):
    document_id = models.IntegerField(blank=True)
    title       = models.TextField(blank=False)
    authors     = models.TextField(blank=False)

    updated     = models.DateTimeField(auto_now=True)
    created     = models.DateTimeField(auto_now_add=True)
    source      = models.CharField(max_length=200, blank=True)

    objects = DocumentManager()

    def __unicode__(self):
        return self.title

    def submitted(self):
        return View.objects.filter(section__document = self).count()

    def available_sections(self):
        return self.section_set.exclude(kind = 'o').all()

    def count_available_sections(self):
        return self.section_set.exclude(kind = 'o').count()

    def get_concepts_for_classification(self):
        # First see if the GM has any annotations for these sections,
        # if not, see if there are any basic concepts for the text
        concepts = ConceptRelationship.objects.filter(annotation__view__user__username = "semmed").values_list('concept', 'target')
        concepts = set(itertools.chain.from_iterable(concepts))
        if len(concepts) >= 2:
          concepts = Concept.objects.in_bulk(concepts).values()
        else:
          concepts = Concept.objects.filter(section__document = doc).distinct()

        concepts = itertools.combinations(concepts, 2)
        concepts = [x for x in concepts]
        random.shuffle(concepts)
        return concepts

    def get_conceptrelation_entries_to_validate(self):
        return ConceptRelationship.objects.filter(
            validate = None,
            annotation__view__section__document = self,
            annotation__view__user__username = "semmed").all()


    def is_complete(self, user, task_type = 'cr'):
        return True if View.objects.filter(user__pk = user.pk, completed = True, task_type = task_type, section__document = self).count() >= self.count_available_sections() else False


    def update_views(self, user, task_type, completed = False):
      for sec in self.available_sections():
          view, created = View.objects.get_or_create(task_type = task_type, section = sec, user = user)
          view.completed = completed
          view.save()


    def annotations(self, username = "goldenmaster", experiment = False):
        if experiment:
          return Annotation.objects.filter(view__section__document = self, view__task_type = "cr", kind = "e", view__user__username=username, experiment = experiment).order_by('start')
        else:
          return Annotation.objects.filter(view__section__document = self, view__task_type = "cr", kind = "e", view__user__username=username).order_by('start')


    def score(self, user):
        gma = self.annotations()

        if user.userprofile.mturk:
          cua = self.annotations(user.username, settings.EXPERIMENT)
        else:
          cua = self.annotations(user.username)


        print " / / / score / / / "
        print gma
        print cua
        compareann = gma[0]
        for ann in gma:
          print ann.is_exact_match(compareann)

        pass

    class Meta:
        ordering = ('created',)



class Section(models.Model):
    SECTION_KIND_CHOICE = (
      ('o', 'Overview'),
      ('t', 'Title'),
      ('a', 'Abstract'),
      ('p', 'Paragraph'),
      ('f', 'Figure'),
    )
    kind = models.CharField(max_length=1, choices=SECTION_KIND_CHOICE)

    text        = models.TextField(blank=True)
    source      = models.ImageField(blank=True, upload_to="media/images/", default = 'images/figure.jpg')

    validate    = models.BooleanField(default = False, blank = True)
    cache       = models.TextField(blank=True)

    concepts    = models.ManyToManyField('Concept', blank=True, null=True)

    updated     = models.DateTimeField(auto_now=True)
    created     = models.DateTimeField(auto_now_add=True)

    document = models.ForeignKey(Document)


    def annotations(self, username = "goldenmaster", experiment = False):
        if experiment:
          return Annotation.objects.filter(view__section = self, view__task_type = "cr", kind = "e", view__user__username=username, experiment = experiment).order_by('start')
        else:
          return Annotation.objects.filter(view__section = self, view__task_type = "cr", kind = "e", view__user__username=username).order_by('start')



    def resultwords(self, user):
        # Gather words and positions from the text
        words_index = WhitespaceTokenizer().span_tokenize(self.text)
        words_text = WhitespaceTokenizer().tokenize(self.text)
        words = zip(words_index, words_text)

        # Add counters for concensus count and personal annotation
        # ((start, stop), 'Word string itself', Intensity, GM Ann ID, User Ann ID, Did user annotate)
        words = [w + (0, None, None, False,) for w in words]

        # Gather other annotations from GM and users for this section
        anns = self.annotations().values_list('pk', 'start', 'text')
        # Build the running counter of times a word was annotated
        for gm_pk, start, text in anns:
          length = len(text)

          for idx, word in enumerate(words):
            word_start = word[0][0]
            counter = word[2]
            if word_start >= start and word_start <= start+length:
              counter += 1
              words[idx] = (word[0], word[1], counter, gm_pk, word[3], word[4])


        if user.userprofile.mturk:
          user_anns = self.annotations(user.username, experiment = settings.EXPERIMENT).values_list('pk', 'start', 'text')
        else:
          user_anns = self.annotations(user.username).values_list('pk', 'start', 'text')
        # Build the running counter of times a word was annotated
        for user_pk, start, text in user_anns:
          length = len(text)
          for idx, word in enumerate(words):
            word_start = word[0][0]
            if word_start >= start and word_start <= start+length:
              words[idx] = (word[0], word[1], word[2], word[3], user_pk, True)

        return words


    def update_view(self, user, task_type, completed = False):
      view = get_object_or_404(View, user = user, task_type = task_type, section = self)
      view.completed = completed
      view.save()
      return view



    def __unicode__(self):
        if self.kind == 'o':
          return '[Overview] '+ self.document.title
        else:
          return self.text


class View(models.Model):
    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    TASK_TYPE_CHOICE = (
      ('cr', 'Concept Recognition'),
      ('cn', 'Concept Normalization'),
      ('rv', 'Relationship Verification'),
      ('ri', 'Relationship Identification'),
      ('rc', 'Relationship Correction'),
    )
    task_type = models.CharField(max_length=3, choices=TASK_TYPE_CHOICE, blank=True, default='cr')
    completed = models.BooleanField(default = False, blank=True)

    section = models.ForeignKey(Section)
    user = models.ForeignKey(User)

    def __unicode__(self):
      return "Doc:"+ str(self.section.document.pk) +", Sec:"+ str(self.section.pk) +" by "+ self.user.username


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
        return "{0} {1}".format(self.message, self.view)


class Comment(models.Model):
    '''
      Refutes without dedicated sections
    '''
    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    message = models.TextField(blank=True)

    TASK_TYPE_CHOICE = (
      ('cr', 'Concept Recognition'),
      ('cn', 'Concept Normalization'),
      ('rv', 'Relationship Verification'),
      ('ri', 'Relationship Identification'),
      ('rc', 'Relationship Correction'),
    )
    task_type = models.CharField(max_length=3, choices=TASK_TYPE_CHOICE, blank=True, default='cr')
    document = models.ForeignKey(Document)
    user = models.ForeignKey(User)

    def __unicode__(self):
        return "{0} {1} {2}".format(self.message, self.document, self.user)




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
    type    = models.CharField(max_length=40, blank=True, null=True)

    text    = models.TextField(blank=True, null=True)
    start   = models.IntegerField(blank=True, null=True)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    user_agent  = models.CharField(max_length=150, blank=True, null=True)
    player_ip   = models.GenericIPAddressField(blank=True, null=True)
    experiment  = models.IntegerField(blank=True, null=True)

    view = models.ForeignKey(View)

    objects = AnnotationManager()


    def simple(self):
      return (self.text, int(self.start))

    def is_exact_match(self, comparing_annotation):
      return True if self.start == comparing_annotation.start and len(self.text) == len(comparing_annotation.text) else False

    def __unicode__(self):
        if self.kind == 'r':
          return "Relationship Ann"
        if self.kind == 'e':
          return "{0} ({1}) [{2}]".format(self.text, self.start, self.pk)
        else:
          return self.type



class Concept(models.Model):
    concept_id = models.TextField(blank=False)
    preferred_name = models.TextField(blank=True)
    definition = models.TextField(blank=True)
    semantic_type = models.TextField(blank=True)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        return "{0} ({1})".format(self.preferred_name, self.concept_id)


class RelationshipType(models.Model):
    full_name = models.CharField(max_length = 80)
    type = models.CharField(max_length = 80)
    definition = models.TextField(blank=True)

    parent = models.ForeignKey("self", blank=True, null=True, related_name="children")

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

        return "{0}:: {1}".format(path, self.full_name)



class ConceptRelationship(models.Model):
    concept       = models.ForeignKey(Concept, related_name="actor")
    relationship  = models.ForeignKey('RelationshipType', blank=True, null=True)
    target        = models.ForeignKey(Concept, blank=True, null=True, related_name="target")

    annotation = models.ForeignKey(Annotation, blank=True, null=True)

    validate = models.ForeignKey("self", blank=True, null=True)
    confidence = models.DecimalField(max_digits=11, decimal_places=5, validators=[MaxValueValidator(1), MinValueValidator(0)], default=0)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        return "{0} >> {1} >> {2}".format(self.concept.preferred_name,
                                          self.relationship.full_name,
                                          self.target)


