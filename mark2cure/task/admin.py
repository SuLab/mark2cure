from django.contrib import admin
from django.db import models

from .models import Task, UserQuestRelationship, DocumentQuestRelationship
from ..common.templatetags.truncatesmart import truncatesmart


class UserQuestRelationshipAdmin(admin.ModelAdmin):
    search_fields = ('user__username', 'task__name')

    list_display = ('task', 'user', 'views_preview',
                    'completed', 'score',
                    'updated', 'created')

    readonly_fields = ('task', 'user', 'views',
                    'completed', 'score',
                    'updated', 'created')

    def views_preview(self, obj):
        return obj.views.count()

    def group(self, obj):
        if obj.task.group:
            return obj.task.group.stub
        else:
            return '[none]'

    mymodel = models.ForeignKey(DocumentQuestRelationship)


class DocumentQuestRelationshipAdmin(admin.ModelAdmin):
    search_fields = ('document__document_id', 'task__name', 'task__group__name')

    list_display = ('task', 'document_preview', 'document_pmid', 'points',
                    'group',
                    'updated', 'created')

    readonly_fields = ('task', 'document',
                       'updated', 'created')

    def document_preview(self, obj):
        return truncatesmart(obj.document.title)

    def document_pmid(self, obj):
        return obj.document.document_id

    def group(self, obj):
        if obj.task.group:
            return obj.task.group.stub
        else:
            return '[none]'

    mymodel = models.ForeignKey(DocumentQuestRelationship)


class TaskAdmin(admin.ModelAdmin):
    search_fields = ('name', 'group__name')

    list_display = ('name', 'contributions', 'completions',
                    'points', 'document_count',
                    'updated', 'created', 'group_preview')

    readonly_fields = ('documents',
                    'users',
                    'updated', 'created', 'group')

    def contributions(self, obj):
        return UserQuestRelationship.objects.filter(task=obj, completed=True).count()

    def group_preview(self, obj):
        if obj.group:
            return obj.group.stub
        else:
            return '[none]'

    def document_count(self, obj):
        return obj.documents.count()

    mymodel = models.ForeignKey(Task)


admin.site.register(Task, TaskAdmin)
admin.site.register(UserQuestRelationship, UserQuestRelationshipAdmin)
admin.site.register(DocumentQuestRelationship, DocumentQuestRelationshipAdmin)
