from django.contrib import admin

from mark2cure.document.models import Document, Pubtator, Section, View, Annotation


class MyInlineModelOptions(admin.TabularInline):
    fields = ('title', 'updated',)
    # define the sortable
    # sortable_field_name = 'position'
    # define sortable_excludes
    sortable_excludes = ('authors',)

admin.site.register(Document)
admin.site.register(Pubtator)
admin.site.register(Section)
admin.site.register(View)
admin.site.register(Annotation)
