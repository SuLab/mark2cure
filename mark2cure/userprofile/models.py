from django.db import models
from django.contrib.auth.models import User

from timezone_field import TimeZoneField
from django_countries.fields import CountryField

from brabeion.models import BadgeAward
from djangoratings.fields import RatingField

import os


class Team(models.Model):
    owner = models.ForeignKey(User)
    name = models.CharField(verbose_name=u'Team Name', help_text=u'You can create a new team.', max_length=255, blank=True)
    created = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        return self.name


def _createHash():
    return os.urandom(40).encode('hex')


def _content_file_name(instance, filename):
    name = _createHash() + os.path.splitext(filename)[1]
    return '/'.join(['avatars', name])


class UserProfile(models.Model):
    user = models.OneToOneField(User, unique=True)
    team = models.ForeignKey(Team, null=True, blank=True)

    timezone = TimeZoneField(default='America/Los_Angeles',
                             blank=True, null=True)
    avatar = models.ImageField(upload_to=_content_file_name,
                               default='images/default.jpg',
                               blank=True)
    rating = RatingField(range=100000, allow_anonymous=True)

    email_notify = models.BooleanField(default=False)

    '''
        Profiling our users
    '''
    MALE = 'm'
    FEMALE = 'f'
    GENDER_CHOICES = (
        (MALE, 'Male'),
        (FEMALE, 'Female'),
    )
    gender = models.CharField(max_length=1,
                              choices=GENDER_CHOICES,
                              blank=True, null=True, default=None)
    age = models.IntegerField(blank=True, null=True, default=None)
    occupation = models.CharField(max_length=255, blank=True)

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
    education = models.IntegerField(choices=EDUCATION_CHOICES,
                                    blank=True, null=True, default=None)
    science_education = models.IntegerField(choices=EDUCATION_CHOICES,
                                            blank=True, null=True,
                                            default=None)
    country = CountryField(blank=True)

    '''
        Profile page features
    '''

    referral = models.TextField(verbose_name=u'I heard about Mark2Cure from', blank=True)
    motivation = models.TextField(verbose_name=u'I contribute to Mark2Cure because', blank=True)
    quote = models.TextField(verbose_name=u'Quote / Signature', blank=True)

    def __unicode__(self):
        return u'Profile of user: %s' % self.user.username

    def highest_level(self, slug='skill'):
        res = BadgeAward.objects.filter(user=self.user, slug=slug).order_by('-level').first()
        return res if res else BadgeAward(name='', level=0)


User.profile = property(lambda u: UserProfile.objects.get_or_create(user=u)[0])
