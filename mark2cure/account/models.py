from django.db import models
from django.db.models import Count

from django.contrib.auth.models import User
from mark2cure.common.models import Message

from mark2cure.document.models import Document, View

from timezone_field import TimeZoneField
import datetime

class UserProfile(models.Model):
    user                  = models.OneToOneField(User, unique=True)

    created_by            = models.ForeignKey(User, null=True, blank=True, related_name="children")
    timezone              = TimeZoneField(default='America/Los_Angeles')

    instructions_enabled  = models.BooleanField(default=True, verbose_name="Display Extra Instructions")

    experience  = models.IntegerField(default=0)
    feedback_0  = models.IntegerField(default=0)
    feedback_1  = models.IntegerField(default=0)
    feedback_2  = models.IntegerField(default=0)
    feedback_3  = models.IntegerField(default=0)

    first_run   = models.BooleanField(default = False, blank = True)
    email_notify  = models.BooleanField(default = False, blank = True)

    mturk           = models.BooleanField(default = False, blank = True)
    turk_last_assignment_id = models.CharField(max_length=200, blank = True)
    turk_submit_to  = models.CharField(max_length=200, blank = True, default = "http://example.com")
    ncbo            = models.BooleanField(default = False, blank = True)


    def score(self, task_type="cr"):
      # Get the last 3 documents with Golden Master Anns
      # documents = Document.objects.annotate(num_annotations=Count('section__view__annotation')).filter(section__view__user__username = "goldenmaster").order_by('-created').distinct()[:3]

      documents = Document.objects.\
        annotate(num_annotations = Count('section__view__annotation')).\
        filter(section__view__user = self, section__view__user__username = "goldenmaster").\
        order_by('-created').\
        distinct()[:3]


      views = View.objects.filter(user = self)


      print " / / / /  DOCS / / / / /"
      print documents
      print " - - "
      for doc in documents:
        print doc.pk, doc.num_annotations

      return documents


    def __unicode__(self):
        return u'Profile of user: %s' % self.user.username


User.profile = property(lambda u: UserProfile.objects.get_or_create(user=u)[0])


class Ncbo(models.Model):
    user = models.OneToOneField(User, unique=True, null=True)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    min_term_size = models.IntegerField()
    score         = models.IntegerField()

    def __unicode__(self):
      return u'NCBO Settings : %s' % self.user.username



User.ncbo = property(lambda u: Ncbo.objects.get_or_create(user=u)[0])
