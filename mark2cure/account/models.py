from django.db import models
from django.db.models import Count
from django.contrib.auth.models import User

from mark2cure.common.models import Message
from mark2cure.document.models import Document, View, Activity

from timezone_field import TimeZoneField
from django_countries.fields import CountryField

import datetime


class UserProfile(models.Model):
    user = models.OneToOneField(User, unique=True)

    created_by = models.ForeignKey(User, null=True, blank = True, related_name = 'children')
    timezone = TimeZoneField(default='America/Los_Angeles')

    instructions_enabled = models.BooleanField(default = True, verbose_name = 'Display Extra Instructions')

    experience = models.IntegerField(default=0)
    feedback_0 = models.IntegerField(default=0)
    feedback_1 = models.IntegerField(default=0)
    feedback_2 = models.IntegerField(default=0)
    feedback_3 = models.IntegerField(default=0)

    first_run = models.BooleanField(default = False, blank = True)
    email_notify = models.BooleanField(default = False, blank = True)

    mturk = models.BooleanField(default = False, blank = True)
    banned = models.BooleanField(default = False, blank = True)
    softblock = models.BooleanField(default = False, blank = True)
    ignore = models.BooleanField(default = False, blank = True)
    current_gm = models.BooleanField(default = False)
    turk_last_assignment_id = models.CharField(max_length=200, blank = True)
    turk_submit_to = models.CharField(max_length = 200, blank = True, default = 'http://example.com')
    ncbo = models.BooleanField(default = False, blank = True)

    user_agent = models.CharField(max_length = 150, blank = True, null = True)
    player_ip = models.GenericIPAddressField(blank = True, null = True)
    '''
        Profiling our users
    '''
    GENDER_CHOICES = (
      ('m', 'Male'),
      ('f', 'Female'),
    )
    gender = models.CharField(max_length = 1, choices = GENDER_CHOICES, blank = True, null = True, default = None)
    age = models.IntegerField(blank = True, null = True, default = None)
    occupation = models.CharField(max_length = 255, blank = True)

    EDUCATION_CHOICES = (
      (0, 'Some elementary'),
      (1, 'Finished elementary'),
      (2, 'Some high school'),
      (3, 'Finished high school'),
      (4, 'Some community college'),
      (5, 'Finished community college'),
      (6, 'Some 4-year college'),
      (7, 'Finished 4-year college'),
      (8, 'Some masters program'),
      (9, 'Finished masters program'),
      (10, 'Some PhD program'),
      (11, 'Finished PhD program'),
    )
    education = models.IntegerField(choices = EDUCATION_CHOICES, blank = True, null = True, default = None)
    science_education = models.IntegerField(choices = EDUCATION_CHOICES, blank = True, null = True, default = None)
    motivation = models.CharField(max_length = 255, blank = True)
    country = CountryField(blank = True)

    def score(self, task_type = 'cr'):
        return sum(Activity.objects.filter(user=self.user, task_type=task_type, submission_type='gm').values_list('f_score', flat=True).all())


    def survey_complete(self):
        if self.gender == None or self.age == None or self.occupation == '' or self.education == None or self.science_education == None or self.motivation == '' or self.country.name == '': return False
        return True


    def __unicode__(self):
        return u'Profile of user: %s' % self.user.username


User.profile = property(lambda u: UserProfile.objects.get_or_create(user=u)[0])


class Ncbo(models.Model):
    user = models.OneToOneField(User, unique = True, null = True)

    updated = models.DateTimeField(auto_now = True)
    created = models.DateTimeField(auto_now_add = True)

    min_term_size = models.IntegerField()
    score = models.IntegerField()

    def __unicode__(self):
      return u'NCBO Settings : %s' % self.user.username


User.ncbo = property(lambda u: Ncbo.objects.get_or_create(user=u)[0])

