from django.contrib import admin
from django.db import models

from mark2cure.common.models import Group, Task, UserQuestRelationship, DocumentQuestRelationship, SupportMessage


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
admin.site.register(DocumentQuestRelationship)
admin.site.register(SupportMessage, SupportMessageAdmin)
