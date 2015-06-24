from mark2cure.common.models import Group, Task, UserQuestRelationship
from mark2cure.userprofile.models import UserProfile, Team

from rest_framework import serializers


class GroupSerializer(serializers.ModelSerializer):
    class Meta:
        model = Group
        fields = ('pk', 'name', 'stub', 'description')


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
        self.user = user
        self.profile = user.profile

        # Instantiate the superclass normally
        super(QuestSerializer, self).__init__(*args, **kwargs)

    user = serializers.SerializerMethodField('get_user_status')
    progress = serializers.SerializerMethodField('get_progress_status')

    def get_user_status(self, task):
        return {'enabled': self.profile.highest_level('skill').level >= task.requires_qualification,
                'completed': UserQuestRelationship.objects.filter(task=task, user=self.user, completed=True).exists()}

    def get_progress_status(self, task):
        current_submissions_count = UserQuestRelationship.objects.filter(task=task, completed=True).count()
        return {'required': task.completions if task.completions else 10000,
                'current': current_submissions_count,
                'completed': current_submissions_count >= task.completions if task.completions else 10000 <= current_submissions_count}

    class Meta:
        model = Task
        fields = ('id', 'name', 'documents', 'points',
                  'requires_qualification', 'provides_qualification',
                  'meta_url', 'user', 'progress')

