from mark2cure.document.models import *
from django.contrib import admin

class MyInlineModelOptions(admin.TabularInline):
    fields = ( "title", "updated",)
    # define the sortable
    # sortable_field_name = "position"
    # define sortable_excludes
    sortable_excludes = ("authors",)

admin.site.register(Document)
admin.site.register(Concept)

admin.site.register(ConceptRelationship)
admin.site.register(RelationshipType)

admin.site.register(Section)
admin.site.register(Refute)
admin.site.register(View)
admin.site.register(Annotation)

