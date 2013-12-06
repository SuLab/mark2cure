from django.db import models
from django.contrib.auth.models import User
from mark2cure.document.models import Document, Annotation

class Concept(models.Model):
    concept_id = models.TextField(unique=True, blank=False)
    annotation = models.ForeignKey(Annotation)


class Message(models.Model):
    message = models.TextField(blank=True)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    user  = models.ForeignKey(User)

    def __unicode__(self):
        return self.message


class Quest(models.Model):
    name    = models.TextField(blank=True)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    documents = models.ManyToManyField(Document)
