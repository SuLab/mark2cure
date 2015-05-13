from django.contrib import admin

from mark2cure.common.models import SkillBadge, PointsBadge, Group, Task, UserQuestRelationship, DocumentQuestRelationship, SupportMessage

admin.site.register(SkillBadge)
admin.site.register(PointsBadge)
admin.site.register(Group)
admin.site.register(Task)
admin.site.register(UserQuestRelationship)
admin.site.register(DocumentQuestRelationship)
admin.site.register(SupportMessage)
