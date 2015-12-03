from mark2cure.userprofile.models import UserProfile, Team
from ..common.models import Group
from ..task.models import Task

from rest_framework import serializers


class GroupSerializer(serializers.ModelSerializer):

    documents = serializers.SerializerMethodField()
    complete_percent = serializers.SerializerMethodField()

    def get_documents(self, group):
        return group.get_documents().values_list('document_id', flat=True)

    def get_complete_percent(self, group):
        return group.percentage_complete()

    class Meta:
        model = Group
        fields = ('pk', 'name', 'stub',
                  'description', 'enabled', 'complete_percent',
                  'documents')


class UserProfileSerializer(serializers.ModelSerializer):

    user = serializers.SerializerMethodField()
    name = serializers.SerializerMethodField()
    score = serializers.SerializerMethodField()

    def get_name(self, profile):
        return profile.user.username

    def get_user(self, profile):
        user = profile.user
        return {'pk': user.pk,
                'username': user.username}

    def get_score(self, profile):
        return int(profile.score) if profile.score else 0

    class Meta:
        model = UserProfile
        fields = ('user', 'name', 'score', 'quote', 'motivation')


class TeamLeaderboardSerializer(serializers.ModelSerializer):

    score = serializers.SerializerMethodField()

    def get_score(self, profile):
        return int(profile.score) if profile.score else 0

    class Meta:
        model = Team
        fields = ('name', 'score',)


class QuestSerializer(serializers.ModelSerializer):

    def __init__(self, *args, **kwargs):
        # Don't pass the 'fields' arg up to the superclass
        context = kwargs.pop('context', {})
        user = context.get('user', None)
        self.user_highest_level = user.profile.highest_level('skill').level

        # Instantiate the superclass normally
        super(QuestSerializer, self).__init__(*args, **kwargs)

    user = serializers.SerializerMethodField('get_user_status')
    progress = serializers.SerializerMethodField('get_progress_status')

    def get_user_status(self, task):
        return {'enabled': self.user_highest_level >= task.requires_qualification,
                'completed': task.user_completed > 0}

    def get_progress_status(self, task):
        return {'required': task.completions if task.completions else 10000,
                'current': task.current_submissions_count,
                'completed': task.current_submissions_count >= task.completions if task.completions else 10000 <= task.current_submissions_count}

    class Meta:
        model = Task
        fields = ('id', 'name', 'documents', 'points',
                  'requires_qualification', 'provides_qualification',
                  'meta_url', 'user', 'progress')

