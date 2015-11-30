from django.contrib import admin
from django.db import models

from .models import Group, SupportMessage
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


class SupportMessageAdmin(admin.ModelAdmin):
    search_fields = ('user__username', 'user__first_name', 'user__last_name',
                     'user__email', 'text')

    list_display = ('get_email', 'text', 'referral',
                    'created')

    readonly_fields = ('user', 'text', 'referral',
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
admin.site.register(SupportMessage, SupportMessageAdmin)
