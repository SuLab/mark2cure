from django.contrib import admin

from ..document.models import Document

from django_comments.admin import CommentsAdmin
from django_comments.models import Comment

from ..common.templatetags.truncatesmart import truncatesmart


class CA(CommentsAdmin):
    list_display = ('text', 'user', 'pmid',
                    'submit_date', 'is_public', 'is_removed')

    readonly_fields = ('content_type', 'object_pk', 'site',
                       'user', 'user_name', 'user_email',
                       'user_url',)

    def text(self, obj):
        return truncatesmart(obj.comment)

    def user(self, obj):
        return obj.user

    def pmid(self, obj):
        return Document.objects.get(pk=obj.object_pk).document_id

admin.site.unregister(Comment)
admin.site.register(Comment, CA)

