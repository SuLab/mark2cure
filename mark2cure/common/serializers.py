from django.contrib.auth.models import User
from mark2cure.common.models import Task, UserQuestRelationship

from rest_framework import serializers


class QuestSerializer(serializers.ModelSerializer):

    def __init__(self, *args, **kwargs):
        # Don't pass the 'fields' arg up to the superclass
        context = kwargs.pop('context', {})
        user = context.get('user', None)
        self.user = user
        self.profile = user.profile

        # Instantiate the superclass normally
        super(QuestSerializer, self).__init__(*args, **kwargs)


    enabled = serializers.SerializerMethodField('get_enabled_status')
    completed = serializers.SerializerMethodField('get_completed_status')
    completions = serializers.SerializerMethodField('get_completion_status')

    def get_enabled_status(self, task):
        return self.profile.highest_level('skill').level >= task.requires_qualification

    def get_completed_status(self, task):
        return UserQuestRelationship.objects.filter(task=task, user=self.user, completed=True).exists()

    def get_completion_status(self, task):
        return 5

    class Meta:
        model = Task
        fields = ('id', 'name', 'documents', 'points',
                  'requires_qualification', 'provides_qualification',
                  'meta_url', 'enabled', 'completed',
                  'completions')

