from django.contrib import admin
from django.db import models

from mark2cure.common.models import Group, Task, UserQuestRelationship, DocumentQuestRelationship, SupportMessage
from mark2cure.common.templatetags.truncatesmart import truncatesmart


class UserQuestRelationshipAdmin(admin.ModelAdmin):

    list_display = ('task', 'user', 'views_preview',
                    'completed', 'score',
                    'updated', 'created')


    readonly_fields = ('task', 'user', 'views',
                    'completed', 'score',
                    'updated', 'created')

    def views_preview(self, obj):
        return obj.views.count()

    def group(self, obj):
        return obj.task.group.stub

    mymodel = models.ForeignKey(DocumentQuestRelationship)


class DocumentQuestRelationshipAdmin(admin.ModelAdmin):

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
        return obj.task.group.stub

    mymodel = models.ForeignKey(DocumentQuestRelationship)


class TaskAdmin(admin.ModelAdmin):

    list_display = ('name', 'kind', 'contributions', 'completions',
                    'points', 'experiment', 'document_count',
                    'requires_qualification', 'provides_qualification', 'meta_url',
                    'updated', 'created', 'group_preview')



    readonly_fields = ('kind', 'documents',
                    'users',
                    'meta_url',
                    'updated', 'created', 'group')

    def contributions(self, obj):
        return UserQuestRelationship.objects.filter(task=obj, completed=True).count()

    def group_preview(self, obj):
        return obj.group.stub

    def document_count(self, obj):
        return obj.documents.count()

    mymodel = models.ForeignKey(Task)


class GroupAdmin(admin.ModelAdmin):

    list_display = ('name', 'enabled', 'stub',
                    'description', 'tasks', 'total_documents')

    def tasks(self, obj):
        return obj.task_set.count()

    def document_count(self, obj):
        return obj.documents.count()

    def total_documents(self, obj):
        return DocumentQuestRelationship.objects.filter(task__group=obj).count()

    mymodel = models.ForeignKey(Task)


class SupportMessageAdmin(admin.ModelAdmin):

    list_display = ( 'get_email', 'text', 'referral',
    'created')

    readonly_fields = ( 'user', 'text', 'referral',
    'created')

    def get_email(self, obj):
        if obj.user:
            return '{first} {last} <{email}>'.format(
                    first=obj.user.first_name,
                    last=obj.user.last_name,
                    email=obj.user.email)
        else:
            return '[anon]'

    mymodel = models.ForeignKey(SupportMessage)


admin.site.register(Group, GroupAdmin)
admin.site.register(Task, TaskAdmin)
admin.site.register(UserQuestRelationship, UserQuestRelationshipAdmin)
admin.site.register(DocumentQuestRelationship, DocumentQuestRelationshipAdmin)
admin.site.register(SupportMessage, SupportMessageAdmin)
