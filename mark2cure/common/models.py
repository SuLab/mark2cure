from django.db import models
from mark2cure.document.models import Document, Annotation

class Concept(models.Model):
  concept_id = models.TextField(unique=True, blank=False)
  annotation = models.ForeignKey(Annotation)

class Message(models.Model):
  message = models.TextField(blank=True)

  updated = models.DateTimeField(auto_now=True)
  created = models.DateTimeField(auto_now_add=True)

class Quest(models.Model):
    name    = models.TextField(blank=True)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    documents = models.ManyToManyField(Document)
