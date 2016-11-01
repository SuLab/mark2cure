from django.contrib import admin
from django.db import models
from django.contrib.humanize.templatetags.humanize import naturaltime

from mark2cure.document.models import Document, Pubtator, Section, View, Annotation
from mark2cure.common.templatetags.truncatesmart import truncatesmart


class DocumentAdmin(admin.ModelAdmin):
    search_fields = ('document_id', 'title', 'source')

    list_display = ('document_id', 'title_preview', 'sections',
            'pubtator', 'annotations', 'completed_views',
            'pending_views', 'source')

    readonly_fields = ('document_id', 'title', 'authors',
            'source')

    def title_preview(self, obj):
        return truncatesmart(obj.title)

    def sections(self, obj):
        return obj.count_available_sections()

    def pubtator(self, obj):
        return obj.valid_pubtator()

    def annotations(self, obj):
        return Annotation.objects.filter(view__section__document=obj).count()

    def completed_views(self, obj):
        return View.objects.filter(section__document=obj, completed=True).count()

    def pending_views(self, obj):
        return View.objects.filter(section__document=obj, completed=False).count()

    mymodel = models.ForeignKey(Document)


class PubtatorAdmin(admin.ModelAdmin):
    search_fields = ('document__document_id', 'session_id', 'content', 'kind')

    list_display = ('pmid', 'kind', 'annotations', 'valid')

    readonly_fields = ('document', 'kind', 'content')

    def pmid(self, obj):
        return obj.document.document_id

    def valid(self, obj):
        return obj.valid() is not False

    def annotations(self, obj):
        return obj.count_annotations()

    mymodel = models.ForeignKey(Pubtator)


class SectionAdmin(admin.ModelAdmin):
    search_fields = ('text', 'document__document_id')

    list_display = ('text_preview', 'kind', 'updated',
            'annotations', 'completed_views', 'pending_views',
            'created', 'document')

    readonly_fields = ('kind', 'text', 'source',
            'updated', 'created', 'document')

    def text_preview(self, obj):
        return truncatesmart(obj.text)

    def annotations(self, obj):
        return Annotation.objects.filter(view__section=obj).count()

    def completed_views(self, obj):
        return View.objects.filter(section=obj, completed=True).count()

    def pending_views(self, obj):
        return View.objects.filter(section=obj, completed=False).count()

    mymodel = models.ForeignKey(Section)


class ViewAdmin(admin.ModelAdmin):
    search_fields = ('user__username', 'opponent__user__username', 'section__document__document_id')

    list_display = ('task_type', 'section_preview', 'annotations',
                    'user', 'partner', 'quest', 'group', 'completed')

    readonly_fields = ('task_type', 'completed', 'opponent',
            'section', 'user',)

    def section_preview(self, obj):
        return truncatesmart(obj.section.text)

    def partner(self, obj):
        if obj.opponent:
            return obj.opponent.user
        else:
            return '[new anns]'

    def annotations(self, obj):
        return Annotation.objects.filter(view=obj).count()

    def quest(self, obj):
        uqr = obj.userquestrelationship_set.first()
        if uqr:
            return uqr.task

    def group(self, obj):
        uqr = obj.userquestrelationship_set.first()
        if uqr:
            if uqr.task:
                if uqr.task.group:
                    return uqr.task.group.stub

    mymodel = models.ForeignKey(View)


class AnnotationAdmin(admin.ModelAdmin):
    search_fields = ('view__user__username', 'view__section__document__document_id')

    list_display = (
            'kind',
            'section', 'username', 'pmid',
            'quest', 'group', 'time_ago',
            'created')

    readonly_fields = ('kind', 'created', 'view')

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
        if uqr:
            if uqr.task:
                return uqr.task

    def group(self, obj):
        uqr = obj.view.userquestrelationship_set.first()
        if uqr:
            if uqr.task:
                if uqr.task.group:
                    return uqr.task.group.stub

    def time_ago(self, obj):
        return naturaltime(obj.created)

    mymodel = models.ForeignKey(Annotation)

admin.site.register(Document, DocumentAdmin)
admin.site.register(Pubtator, PubtatorAdmin)
admin.site.register(Section, SectionAdmin)
admin.site.register(View, ViewAdmin)
admin.site.register(Annotation, AnnotationAdmin)
