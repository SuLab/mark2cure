from django.contrib.auth.models import User
from mark2cure.common.models import Task, UserQuestRelationship

from rest_framework import serializers
from random import randrange

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
        return {'required': task.completions,
                'current': current_submissions_count,
                'completed': task.completions == current_submissions_count}

    class Meta:
        model = Task
        fields = ('id', 'name', 'documents', 'points',
                  'requires_qualification', 'provides_qualification',
                  'meta_url', 'user', 'progress')

