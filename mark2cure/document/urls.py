from django.conf.urls import patterns, include, url

from rest_framework import routers

from mark2cure.document.views import *


router = routers.DefaultRouter()

urlpatterns = patterns(
    'mark2cure.document.views',

    url(r'^(?P<task_id>\d+)/(?P<doc_id>\d+)/$', r'identify_annotations'),
    url(r'^(?P<task_id>\d+)/(?P<doc_id>\d+)/results/$', r'identify_annotations_results'),
    url(r'^(?P<task_id>\d+)/(?P<doc_id>\d+)/section/(?P<section_id>\d+)/annotation/create/$',
        r'identify_annotations_submit'),

    url(r'^(?P<task_id>\d+)(?P<doc_id>\d+)/submit/$', r'submit'),

    # REST Framework
    url(r'^(?P<doc_id>\d+)/section/(?P<section_id>\d+)/results/top/$',
        TopUserViewSet.as_view()),
    url('^(?P<doc_id>\d+)/section/(?P<section_id>\d+)/user/(?P<user_id>\d+)/annotations/$',
        AnnotationViewSet.as_view()),
    url(r'^', include(router.urls)),
)
