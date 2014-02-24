from django.db import models
from django.conf import settings
from django.core.validators import MaxValueValidator, MinValueValidator
from django.utils.encoding import smart_text

from mark2cure.document.managers import DocumentManager, AnnotationManager
from django.contrib.auth.models import User

from ttp import ttp
from decimal import Decimal as D
from copy import copy

import requests, random, datetime

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

    def __unicode__(self):
        if self.kind == 'o':
          return '[Overview] '+ self.document.title
        else:
          return self.text


class View(models.Model):
    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    section = models.ForeignKey(Section)
    user = models.ForeignKey(User)

    def __unicode__(self):
      return "Doc:"+ str(self.section.document.pk) +", Sec:"+ str(self.section.pk) +" by "+ self.user.username


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

    def __unicode__(self):
        if self.kind == 'r':
          return "Relationship Ann"
        else:
          return self.type

    def simple(self):
      return (self.text, int(self.start))



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

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        return "{0} >> {1} >> {2}".format(self.concept.preferred_name,
                                          self.relationship.full_name,
                                          self.target)



