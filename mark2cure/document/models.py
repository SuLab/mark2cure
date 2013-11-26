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

class Annotation(models.Model):
    ANNOTATION_KIND_CHOICE = (
      ('e', 'Entities'),
      ('a', 'Attributes'),
      ('r', 'Relations'),
      ('t', 'Triggers'),
      ('e', 'Events'),
    )
    kind = models.CharField(max_length=1, choices=ANNOTATION_KIND_CHOICE)

    # Disease, Gene, Protein, et cetera...
    type    = models.CharField(max_length=12, blank=True)

    text    = models.TextField(blank=False)
    start   = models.IntegerField()

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    user_agent  = models.CharField(max_length=150, blank=True)
    player_ip   = models.GenericIPAddressField()
    experiment  = models.IntegerField()


class View(models.Model):
    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    user = models.ForeignKey(User)
    annotation = models.ForeignKey(Annotation)


class Section(models.Model):
    SECTION_KIND_CHOICE = (
      ('p', 'Paragraph'),
      ('f', 'Figure'),
      ('t', 'Table'),
    )
    kind = models.CharField(max_length=1, choices=SECTION_KIND_CHOICE)

    text        = models.TextField(blank=False)
    source      = models.ImageField(blank=True, upload_to="media/images/", default = 'images/figure.jpg')

    validate    = models.BooleanField(default = False, blank = True)
    cache       = models.TextField(blank=True)

    view = models.ForeignKey(View)


class Document(models.Model):
    document_id = models.IntegerField(blank=True)
    title       = models.TextField(unique=True, blank=False)
    authors     = models.TextField(blank=False)

    updated     = models.DateTimeField(auto_now=True)
    created     = models.DateTimeField(auto_now_add=True)
    source      = models.CharField(max_length=200, blank=True)

    section = models.ForeignKey(Section)

    objects = DocumentManager()

    def __unicode__(self):
        return self.title


# def json_view(self, user):
#     if user.is_anonymous():
#       viewed = False
#       annotations = []
#     else:
#       viewed = True if db.session.query(View).filter_by(user = user).filter_by(document = self).first() else False
#       annotations = db.session.query(Annotation).filter_by(user = user).filter_by(document = self).all()

#     popularity = self.cache.split(", ") if self.cache else ["0"]

#     return {  'id'          : self.id,
#               'document_id' : self.document_id,
#               'text'        : self.text,
#               'title'       : self.title,
#               'created'     : self.created.isoformat(),
#               # relationships
#               'complete'    : viewed,
#               'annotations' : [i.json_view() for i in annotations],
#               'popularity'  : popularity
#               }

# def compare_view(self):
#   # Returns back the text dictionary for comparision
#   offset = len(self.text) - len(self.text.lstrip())
#   return {
#             'text'      : self.text.strip().lower(),
#             'start'     : int(self.start)+offset,
#           }
