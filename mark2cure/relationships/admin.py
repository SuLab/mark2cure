'''
Jennifer new admin page (involved in registering changes to the database)
'''
from django.contrib import admin
from django.db import models


from mark2cure.relationships.models import Paper

# TODO
'''
look into search_fields, list_display, readonly_fields
'''

# Overall abstract text basics ()
"""
class PaperAdmin(admin.ModelAdmin):
    fields = ['pmid', 'title', 'abstract', 'annotations']

"""


class PaperAdmin(admin.ModelAdmin):
    #list_display = ('pmid', 'title', 'abstract', 'annotations', 'chemicals', 'diseases', 'gold_relations')
    list_display = ('pmid', 'title', 'abstract', 'annotations', 'relations')

    # TODO add defs here similar to max's admin models


admin.site.register(Paper, PaperAdmin)


"""
    # TODO is there a title for the documents Toby uses?
    list_display = ('document_id', 'pmid', 'title_text', 'abstract_text')

# sentences (should not be stored as entire sentences), but should be stored as
# ingegers with information about the "span" of each sentence.

"""
"""
class SectionAdmin(admin.ModelAdmin):


    # TODO document should be the ForeignKey
    list_display = ('sentence_span', 'chemicals', 'diseases')



    # logic: sentence span, chemical (c) span, disease (d) span.
    # Does c or d fall into the same span as the sentence?
    # unique combination of chemical/disease related to the sentence, abstract,

class RelationAdmin(admin.ModelAdmin):
"""
