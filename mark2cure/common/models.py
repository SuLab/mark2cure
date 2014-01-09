from django.db import models
from django.contrib.auth.models import User

class Concept(models.Model):
    concept_id = models.TextField(blank=False)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        return self.concept_id


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
