from ..userprofile.models import Team
from ..common.models import Group
from ..task.models import Task

from rest_framework import serializers
from django.contrib.auth.models import User


class NERGroupSerializer(serializers.ModelSerializer):

    complete_percent = serializers.SerializerMethodField()

    def get_complete_percent(self, group):
        return group.percentage_complete()

    class Meta:
        model = Group
        fields = ('pk', 'name', 'stub',
                  'description', 'enabled', 'complete_percent')


class LeaderboardSerializer(serializers.ModelSerializer):

    user = serializers.SerializerMethodField()
    name = serializers.SerializerMethodField()
    score = serializers.SerializerMethodField()

    def get_name(self, user):
        return user.username

    def get_user(self, user):
        return {'pk': user.pk,
                'username': user.username}

    def get_score(self, user):
        return int(user.score) if user.score else 0

    class Meta:
        model = User
        fields = ('user', 'name', 'score',)


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
        # context = kwargs.pop('context', {})
        # user = context.get('user', None)

        # Instantiate the superclass normally
        super(QuestSerializer, self).__init__(*args, **kwargs)

    user = serializers.SerializerMethodField('get_user_status')
    progress = serializers.SerializerMethodField('get_progress_status')

    def get_user_status(self, task):
        if hasattr(task, 'user_completed'):
            # All NER Tasks require the same qualiciation: has or has not completed all NER Training
            # We used to segment this further but is no longer used so leaving in for historic purposes
            # until positive all references to it are removed
            return {'enabled': True,
                    'completed': task.user_completed > 0}
        else:
            return {'enabled': False,
                    'completed': False}

    def get_progress_status(self, task):
        return {'required': task.completions if task.completions else 10000,
                'current': task.current_submissions_count,
                'completed': task.current_submissions_count >= task.completions if task.completions else 10000 <= task.current_submissions_count}

    class Meta:
        model = Task
        fields = ('id', 'name', 'documents', 'points',
                  'user', 'progress')


class DocumentRelationSerializer(serializers.Serializer):

    id = serializers.IntegerField()
    document_id = serializers.IntegerField()
    title = serializers.CharField()

    total_document_relationships = serializers.IntegerField()
    user_document_relationships = serializers.IntegerField()

    community_completed = serializers.BooleanField()
    community_progress = serializers.FloatField()

    user_completed = serializers.BooleanField()
    user_progress = serializers.FloatField()
    user_answered = serializers.IntegerField()

