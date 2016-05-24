from django.db import models
from django.contrib.auth.models import User
from django.contrib.contenttypes.models import ContentType
from django.conf import settings
from django.db.models import Q

from timezone_field import TimeZoneField
from django_countries.fields import CountryField

from djangoratings.fields import RatingField

from ..document.models import Annotation, Document, View
from ..common.models import Group
from ..task.models import Task, UserQuestRelationship, Level
from ..task.relation.models import RelationAnnotation
from ..task.entity_recognition.models import EntityRecognitionAnnotation
from ..analysis.models import Report
from ..score.models import Point

from django.utils import timezone
import pandas as pd
import datetime
import os


class Team(models.Model):
    owner = models.ForeignKey(User)
    name = models.CharField(verbose_name=u'Team Name', help_text=u'You can create a new team.', max_length=255, blank=True)
    description = models.TextField(blank=True, null=True)
    created = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        return self.name

    def last_active_user(self):
        member_profile = self.userprofile_set.order_by('-last_seen').first()
        if member_profile:
            return member_profile.user
        else:
            return None

    def members_count(self):
        return self.userprofile_set.count()

    def total_annotations(self):
        team_user_ids = self.userprofile_set.values_list('user__pk', flat=True)
        return Annotation.objects.filter(view__user__pk__in=team_user_ids)

    def total_documents(self):
        doc_pks = set(self.total_annotations().values_list('view__section__document__pk', flat=True))
        return Document.objects.filter(pk__in=doc_pks)

    def finished_quests(self):
        team_user_ids = self.userprofile_set.values_list('user__pk', flat=True)
        return UserQuestRelationship.objects.filter(user__pk__in=team_user_ids, completed=True)

    def total_score(self):
        # (TODO) user.id vs userprofile.id checks here
        from mark2cure.api.views import users_with_score
        users_queryset = users_with_score(days=10000)
        team_user_pks = self.userprofile_set.values_list('pk', flat=True)
        return sum(filter(None, [row.score for row in filter(lambda x: x.id in team_user_pks, users_queryset)]))

    def current_avg_f(self, weighted=True):
        '''
            Return back the weighted mean (pairings count) f-score
        '''
        reports = Report.objects.filter(report_type=1).order_by('-created')[:Group.objects.count()]

        team_user_profile_pks = self.userprofile_set.values_list('pk', flat=True)
        team_user_profile_pks = [str(u) for u in team_user_profile_pks]

        dataframes = [r.dataframe[r.dataframe['user'].isin(team_user_profile_pks)] for r in reports]

        try:
            team_df = pd.concat(dataframes)

            if weighted:
                team_df['wf'] = team_df['pairings'] * team_df['f-score']
                return team_df['wf'].sum() / team_df['pairings'].sum()
            else:
                return team_df['f-score'].mean()
        except:
            return 0.0


def _createHash():
    return os.urandom(40).encode('hex')


def _content_file_name(instance, filename):
    name = _createHash() + os.path.splitext(filename)[1]
    return '/'.join(['avatars', name])


class UserProfile(models.Model):
    user = models.OneToOneField(User, unique=True)
    team = models.ForeignKey(Team, null=True, blank=True)
    last_seen = models.DateTimeField(null=True, blank=True)

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

    def score(self, task=None):

        if task:
            if 'entity' in task:
                val = sum(Point.objects.filter(user=self.user).filter(
                    Q(object_id__isnull=True) | Q(content_type=ContentType.objects.get_for_model(EntityRecognitionAnnotation))
                ).values_list('amount', flat=True))
                print val

                # These are being assigned to Tasks
                # (TODO) This whole thing is fucked up, need to get away from Tasks and
                # find clear way to allocate units complete and points for multi-game system
                quest_task_pks = Task.objects.filter(kind='q').values_list('pk', flat=True)
                val += sum(Point.objects.filter(
                    user=self.user,
                    content_type=ContentType.objects.get_for_model(Task),
                    object_id__in=quest_task_pks
                ).values_list('amount', flat=True))

            if 'relation' in task:

                # Points for submitting individual relation steps
                val = sum(Point.objects.filter(
                    user=self.user,
                    content_type=ContentType.objects.get_for_model(RelationAnnotation)
                ).values_list('amount', flat=True))

                # Points for submitting relation sets
                completed_relation_view_pks = View.objects.filter(user=self.user, task_type='ri').values_list('pk', flat=True)
                val += sum(Point.objects.filter(
                    user=self.user,
                    content_type=ContentType.objects.get_for_model(View),
                    object_id__in=completed_relation_view_pks
                ).values_list('amount', flat=True))

        else:
            val = sum(Point.objects.filter(user=self.user).values_list('amount', flat=True))

        return int(val)

    def unlocked_tasks(self):
        arr = []
        if Level.objects.filter(user=self.user, task_type='e', level=7).exists():
            arr.append('entity_recognition')

        if Level.objects.filter(user=self.user, task_type='r', level=3).exists():
            arr.append('relation')

        return arr

    def all_annotations(self):
        return Annotation.objects.filter(view__user=self.user)

    def quests_count(self):
        return UserQuestRelationship.objects.filter(user=self.user, completed=True).count()

    def completed_document_pks(self):
        return list(set(View.objects.filter(user=self.user, completed=True).values_list('section__document', flat=True)))

    def contributed_groups(self):
        views = View.objects.filter(
            user=self.user,
            completed=True).values_list('userquestrelationship__task__group', flat=True)
        group_pks = list(set(views))
        group_pks = [x for x in group_pks if x is not None]
        return Group.objects.filter(pk__in=group_pks).all()

    def current_avg_f(self, weighted=True):
        '''
            Return back the weighted mean (pairings count) f-score
        '''
        reports = Report.objects.filter(report_type=1).order_by('-created')[:Group.objects.count()]
        dataframes = [r.dataframe[r.dataframe['user'] == str(self.user.pk)] for r in reports]

        try:
            user_df = pd.concat(dataframes)
            if weighted:
                user_df['wf'] = user_df['pairings'] * user_df['f-score']
                return user_df['wf'].sum() / user_df['pairings'].sum()
            else:
                return user_df['f-score'].mean()
        except:
            return 0.0

    def online(self):
        if self.last_seen:
            if timezone.now() > self.last_seen + datetime.timedelta(seconds=settings.USER_ONLINE_TIMEOUT):
                return False
            else:
                return True
        else:
            return False

    def available_quests(self):
        '''
            This returns back the contept of how many quests are available for
            the specific user.

            It repects community completion in addition to user completion

            Returns: Units of Quests available

        '''

        # Quests which are available and uncompleted
        # by the user, annotated by how many times the community
        # has completed them
        queryset = Task.objects.filter(
            kind=Task.QUEST,
            group__enabled=True,
        ).extra(select={
            "community_submission_count": """
                SELECT COUNT(*) AS cummunity_submission_count
                FROM task_userquestrelationship
                WHERE (task_userquestrelationship.completed = 1
                    AND task_userquestrelationship.task_id = task_task.id)""",

            "user_completed": """
                SELECT COUNT(*) AS user_completed
                FROM task_userquestrelationship
                WHERE (task_userquestrelationship.completed = 1
                    AND task_userquestrelationship.user_id = %d
                    AND task_userquestrelationship.task_id = task_task.id)""" % (self.user.pk,)

        }).values('id', 'community_submission_count', 'user_completed', 'completions',)
        uncompleted_quests = [task for task in queryset if task['user_completed'] == 0 and (task['completions'] is None or task['community_submission_count'] < task['completions'])]
        return len(uncompleted_quests)


User.profile = property(lambda u: UserProfile.objects.get_or_create(user=u)[0])
