from mark2cure.document.models import Document, Section, View, Annotation
from django.contrib import admin

class MyInlineModelOptions(admin.TabularInline):
    fields = ( "title", "updated",)
    # define the sortable
    # sortable_field_name = "position"
    # define sortable_excludes
    sortable_excludes = ("authors",)

admin.site.register(Document)
admin.site.register(Section)
admin.site.register(View)
admin.site.register(Annotation)


