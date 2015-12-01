
from django.contrib import admin


from mark2cure.relation.models import Concept, Relation, Answer

class ConceptAdmin(admin.ModelAdmin):
    list_display = ('document', 'uid', 'stype', 'text')

class RelationAdmin(admin.ModelAdmin):
    list_display = ('document', 'relation', 'concept1_id', 'concept2_id')

class AnswerAdmin(admin.ModelAdmin):
    list_display = ('relation', 'relation_type', 'username')

admin.site.register(Concept, ConceptAdmin)
admin.site.register(Relation, RelationAdmin)
admin.site.register(Answer, AnswerAdmin)
