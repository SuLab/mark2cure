'''
doc namespaced URLs
'''

from django.conf import settings
from django.conf.urls import patterns, include, url

urlpatterns = patterns('mark2cure.document.views',
    url(r'^$', r'list'),
    url(r'^page/(?P<page_num>\d+)/$', r'list'),
    url(r'^(?P<doc_id>\d+)/$', r'read'),
    url(r'^(?P<doc_id>\d+)/annotations/(?P<section_id>\d+)$', r'submit_annotations'),
    url(r'^create/$', r'create'),
)
