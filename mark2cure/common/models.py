from django.db import models
from django.contrib.auth.models import User


class Message(models.Model):
    message = models.TextField(blank=True)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    user  = models.ForeignKey(User)

    def __unicode__(self):
        return self.message


class SurveyFeedback(models.Model):

    question = models.CharField(max_length=40, blank=True, null=True)
    response = models.TextField(max_length=40, blank=True, null=True)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    user  = models.ForeignKey(User)

    def __unicode__(self):
      return "{0}: {1} ({2})".format(self.question, self.response, self.user.username)


class Quest(models.Model):
    name    = models.TextField(blank=True)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)


