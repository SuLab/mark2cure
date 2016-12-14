from django.conf import settings

from ..userprofile.models import Team
from ..common.models import Document, Group
from ..task.models import Level, Task

from rest_framework import serializers
from django.contrib.auth.models import User


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
        context = kwargs.pop('context', {})
        user = context.get('user', None)
        if user.is_authenticated() and Level.objects.filter(user=user, task_type='e').exists():
            self.user_highest_level = Level.objects.filter(user=user, task_type='e').first().level
        else:
            self.user_highest_level = 0

        # Instantiate the superclass normally
        super(QuestSerializer, self).__init__(*args, **kwargs)

    user = serializers.SerializerMethodField('get_user_status')
    progress = serializers.SerializerMethodField('get_progress_status')

    def get_user_status(self, task):
        if hasattr(task, 'user_completed'):
            return {'enabled': self.user_highest_level >= task.requires_qualification,
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
                  'requires_qualification', 'provides_qualification',
                  'meta_url', 'user', 'progress')


class DocumentRelationSerializer(serializers.Serializer):

    id = serializers.IntegerField()
    document_id = serializers.IntegerField()
    title = serializers.CharField()
    relationships = serializers.IntegerField()
    progress = serializers.FloatField()

