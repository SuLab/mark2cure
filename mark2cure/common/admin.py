from django.contrib import admin
from django.db import models

from .models import Group
from ..task.models import DocumentQuestRelationship


class GroupAdmin(admin.ModelAdmin):
    search_fields = ['name', 'description', 'stub']

    list_display = ('name', 'order', 'mean_f',
                    'enabled', 'stub', 'description',
                    'tasks', 'total_documents')

    def tasks(self, obj):
        return obj.task_set.count()

    def mean_f(self, obj):
        return obj.current_avg_f()

    def total_documents(self, obj):
        return DocumentQuestRelationship.objects.filter(task__group=obj).count()

    mymodel = models.ForeignKey(Group)


admin.site.register(Group, GroupAdmin)
