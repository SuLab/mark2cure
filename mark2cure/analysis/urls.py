'''
analysis namespaced URLs
'''

from django.conf import settings
from django.conf.urls import patterns, include, url
from rest_framework import routers

from mark2cure.document.views import *

router = routers.DefaultRouter()
router.register(r'relationshiptypes', RelationshipTypeViewSet)

urlpatterns = patterns('mark2cure.analysis.views',
    url(r'^network$', r'network'),
    # url(r'^page/(?P<page_num>\d+)/$', r'list'),
)

# //-- Specific Doc for relationships
# 'document/:doc_id/relationship' : 'showDocRelationship',
#
# //-- Specific Document
# 'document/:doc_id?assignmentId=:var1&hitId=:var2&workerId=:var3&turkSubmitTo=:var4' : 'showDocument',
# 'document/:doc_id?assignmentId=:var1&hitId=:var2' : 'showDocument',
# 'document/:doc_id' : 'showDocument',

