'''
doc namespaced URLs
'''

from django.conf import settings
from django.conf.urls import patterns, include, url
from rest_framework import routers

from mark2cure.document.views import *

router = routers.DefaultRouter()
router.register(r'relationshiptypes', RelationshipTypeViewSet)

urlpatterns = patterns('mark2cure.document.views',
    url(r'^$', r'list'),
    url(r'^page/(?P<page_num>\d+)/$', r'list'),

    url(r'^(?P<doc_id>\d+)/$', r'identify_annotations'),
    url(r'^(?P<doc_id>\d+)/results/$', r'identify_annotations_results'),
    url(r'^(?P<doc_id>\d+)/section/(?P<section_id>\d+)/annotation/create/$', r'identify_annotations_submit'),
    url(r'^(?P<doc_id>\d+)/section/(?P<section_id>\d+)/refute/create/$', r'refute_section'),

    url(r'^(?P<doc_id>\d+)/concepts/validate/$', r'validate_concepts'),
    url(r'^(?P<doc_id>\d+)/concepts/validate/submit/$', r'validate_concepts_submit'),

    url(r'^(?P<doc_id>\d+)/concepts/identify/$', r'identify_concepts'),
    url(r'^(?P<doc_id>\d+)/concepts/identify/submit/$', r'identify_concepts_submit'),

    url(r'^(?P<doc_id>\d+)/submit/$', r'submit'),
    url(r'^(?P<doc_id>\d+)/next/$', r'next'),
    url(r'^create/$', r'create'),
    url(r'^(?P<doc_id>\d+)/delete/$', r'create'),

    # REST Framework
    url(r'^', include(router.urls)),

)

# //-- Specific Doc for relationships
# 'document/:doc_id/relationship' : 'showDocRelationship',
#
# //-- Specific Document
# 'document/:doc_id?assignmentId=:var1&hitId=:var2&workerId=:var3&turkSubmitTo=:var4' : 'showDocument',
# 'document/:doc_id?assignmentId=:var1&hitId=:var2' : 'showDocument',
# 'document/:doc_id' : 'showDocument',
