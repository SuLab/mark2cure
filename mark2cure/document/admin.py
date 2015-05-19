from django.contrib import admin
from django.db import models
from django.contrib.humanize.templatetags.humanize import naturaltime

from mark2cure.document.models import Document, Pubtator, Section, View, Annotation


class MyInlineModelOptions(admin.TabularInline):
    fields = ('title', 'updated',)
    # define the sortable
    # sortable_field_name = 'position'
    # define sortable_excludes
    sortable_excludes = ('authors',)

class AnnotationAdmin(admin.ModelAdmin):
    list_display = ('text', 'type', 'start',
            'section', 'username', 'pmid',
            'quest', 'group', 'time_ago',
            'created')

    readonly_fields = ('kind', 'type', 'text',
            'start', 'created', 'view')

    def section(self, obj):
        if obj.view:
            return obj.view.section.get_kind_display()

    def username(self, obj):
        if obj.view.user:
            return obj.view.user.username

    def pmid(self, obj):
        if obj.view.section:
            return obj.view.section.document.document_id

    def quest(self, obj):
        uqr = obj.view.userquestrelationship_set.first()
        return uqr.task

    def group(self, obj):
        uqr = obj.view.userquestrelationship_set.first()
        return uqr.task.group.stub

    def time_ago(self, obj):
        return naturaltime(obj.created)

    mymodel = models.ForeignKey(Annotation)

admin.site.register(Document)
admin.site.register(Pubtator)
admin.site.register(Section)
admin.site.register(View)
admin.site.register(Annotation, AnnotationAdmin)
