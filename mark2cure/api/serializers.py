from mark2cure.userprofile.models import UserProfile, Team
from ..common.models import Document, Group
from ..task.models import Task

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


class UserSerializer(serializers.ModelSerializer):

    user = serializers.SerializerMethodField()
    name = serializers.SerializerMethodField()
    score = serializers.SerializerMethodField()
    quote = serializers.SerializerMethodField()
    motivation = serializers.SerializerMethodField()

    def get_name(self, user):
        return user.username

    def get_user(self, user):
        return {'pk': user.pk,
                'username': user.username}

    def get_score(self, user):
        return int(user.score) if user.score else 0

    def get_quote(self, user):
        return user.profile.quote

    def get_motivation(self, user):
        return user.profile.motivation

    class Meta:
        model = User
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


class DocumentRelationSerializer(serializers.ModelSerializer):

    def __init__(self, *args, **kwargs):
        # Don't pass the 'fields' arg up to the superclass
        context = kwargs.pop('context', {})
        user = context.get('user', None)

        # Instantiate the superclass normally
        super(DocumentRelationSerializer, self).__init__(*args, **kwargs)

    user = serializers.SerializerMethodField('get_user_status')
    progress = serializers.SerializerMethodField('get_progress_status')

    def get_user_status(self, document):
        return {'enabled': True,
                'completed': True if document.user_completed_count > 0 else False}

    def get_progress_status(self, task):
        return {'required': 15,
                'current': task.current_completed_count,
                'completed': task.current_completed_count >= 15 }

    class Meta:
        model = Document
        fields = ('id', 'title', 'document_id',
                  'user', 'progress')

