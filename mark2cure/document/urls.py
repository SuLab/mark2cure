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
    url(r'^(?P<doc_id>\d+)/comment/create/$', r'comment_document'),

    url(r'^(?P<doc_id>\d+)/submit/$', r'submit'),
    url(r'^(?P<doc_id>\d+)/next/$', r'next'),
    url(r'^create/$', r'create'),
    url(r'^(?P<doc_id>\d+)/delete/$', r'create'),

    # REST Framework
    url(r'^(?P<doc_id>\d+)/section/(?P<section_id>\d+)/results/top/$', TopUserViewSet.as_view()),
    url('^(?P<doc_id>\d+)/section/(?P<section_id>\d+)/user/(?P<user_id>\d+)/annotations/$', AnnotationViewSet.as_view()),
    url(r'^', include(router.urls)),
    )

