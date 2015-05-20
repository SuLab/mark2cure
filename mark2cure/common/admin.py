from django.contrib import admin
from django.db import models

from mark2cure.common.models import Group, Task, UserQuestRelationship, DocumentQuestRelationship, SupportMessage
from mark2cure.common.templatetags.truncatesmart import truncatesmart


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


admin.site.register(Group)
admin.site.register(Task)
admin.site.register(UserQuestRelationship)
admin.site.register(DocumentQuestRelationship, DocumentQuestRelationshipAdmin)
admin.site.register(SupportMessage, SupportMessageAdmin)
