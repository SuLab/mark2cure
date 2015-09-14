'''
Jennifer new admin page (involved in registering changes to the database)
'''
from django.contrib import admin


from mark2cure.relation.models import Paper, Annotation, Sentence, Relation, Answer

# TODO
'''
look into search_fields, list_display, readonly_fields
'''


class PaperAdmin(admin.ModelAdmin):
    list_display = ('pmid', 'title', 'abstract')

    # TODO add defs here similar to max's admin models? or in tasks?


class AnnotationAdmin(admin.ModelAdmin):
    list_display = ('document', 'uid', 'stype', 'text', 'start', 'stop')


class SentenceAdmin(admin.ModelAdmin):
    list_display = ('paper', 'uid', 'text', 'start', 'stop', 'annotations')


class RelationAdmin(admin.ModelAdmin):
    list_display = ('document', 'relation', 'chemical_id', 'disease_id', 'automated_cid')


class AnswerAdmin(admin.ModelAdmin):
    list_display = ('relation', 'relation_pair', 'relation_type', 'user_confidence', 'username')


admin.site.register(Paper, PaperAdmin)
admin.site.register(Annotation, AnnotationAdmin)
admin.site.register(Sentence, SentenceAdmin)
admin.site.register(Relation, RelationAdmin)
admin.site.register(Answer, AnswerAdmin)

"""
# sentences (should not be stored as entire sentences), but should be stored as
# ingegers with information about the "span" of each sentence.

"""
