from django.db import models
from django.conf import settings
from django.core.validators import MaxValueValidator, MinValueValidator
from django.utils.encoding import smart_text

from mark2cure.document.managers import DocumentManager, AnnotationManager

from ttp import ttp
from decimal import Decimal as D
from copy import copy

import requests, random, datetime

class Document(models.Model):
    document_id = models.IntegerField()
    text        = models.TextField(blank=False)
    title       = models.TextField(blank=False)

    updated     = models.DateTimeField(auto_now=True)
    created     = models.DateTimeField(auto_now_add=True)

    cache       = models.TextField(blank=True)
    source      = models.CharField(max_length=200, blank=True)

    validate    = models.BooleanField(default = False, blank = True)

    # views           = db.relationship('View',             backref=db.backref('document',  lazy='select'))
    # annotations     = db.relationship('Annotation',       backref=db.backref('document',  lazy='select'))
    # quest_relations = db.relationship('QuestRelation', backref=db.backref('document',  lazy='select'))


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




class Annotation(models.Model):
    kind    = models.IntegerField()
    type    = models.CharField(max_length=12, blank=True)
    text    = models.TextField(blank=False)
    start   = models.IntegerField()
    length  = models.IntegerField()
    stop    = models.IntegerField()

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    user_agent  = models.CharField(max_length=150, blank=True)
    player_ip   = models.GenericIPAddressField()
    experiment  = models.IntegerField()

    # user_id       = db.Column(db.Integer, db.ForeignKey('user.id'))
    # document_id   = db.Column(db.Integer, db.ForeignKey('document.id'))
    # concept_id    = db.Column(db.Integer, db.ForeignKey('concept.id'))

    def __unicode__(self):
        return self.text


    def compare_view(self):
      # Returns back the text dictionary for comparision
      offset = len(self.text) - len(self.text.lstrip())
      return {
                'text'      : self.text.strip().lower(),
                'start'     : int(self.start)+offset,
              }


class View(models.Model):
    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    # user_id       = db.Column(db.Integer, db.ForeignKey('user.id'))
    # document_id   = db.Column(db.Integer, db.ForeignKey('document.id'))
