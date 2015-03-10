from django.conf.urls import patterns, include, url

from rest_framework import routers

from . import views

router = routers.DefaultRouter()

urlpatterns = patterns('',

    url(r'^(?P<pubmed_id>\d+).(?P<format>\w+)$',
        views.read_pubmed_bioc, name='read-pubmed-bioc'),

    url(r'^(?P<task_id>\d+)/(?P<doc_id>\d+)/$',
        views.identify_annotations, name='read'),

    url(r'^(?P<task_id>\d+)/(?P<doc_id>\d+)/results/$',
        views.identify_annotations_results, name='results'),

    url(r'^(?P<task_id>\d+)/(?P<doc_id>\d+)/section/(?P<section_id>\d+)/annotation/create/$',
        views.identify_annotations_submit, name='create'),

    url(r'^(?P<task_id>\d+)/(?P<doc_id>\d+)/submit/$',
        views.submit, name='submit'),

    # REST Framework
    url('^(?P<doc_id>\d+)/section/(?P<section_id>\d+)/user/(?P<user_id>\d+)/annotations/$',
        views.AnnotationViewSet.as_view()),
    url(r'^', include(router.urls)),
)
