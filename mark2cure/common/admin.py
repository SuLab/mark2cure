from django.contrib import admin

from mark2cure.common.models import Group, Task, UserQuestRelationship, DocumentQuestRelationship, SupportMessage

admin.site.register(Group)
admin.site.register(Task)
admin.site.register(UserQuestRelationship)
admin.site.register(DocumentQuestRelationship)
admin.site.register(SupportMessage)
