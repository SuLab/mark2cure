from django.conf.urls import patterns, url

from . import views

urlpatterns = patterns('',


    url(r'^(?P<document_pk>\d+)/(?P<relation_pk>\d+)/submit/$',
        views.submit_annotation, name='submit-annotation'),

    url(r'^(?P<document_pk>\d+)/results/$',
        views.show_document_results, name='results'),

    url(r'^(?P<document_pk>\d+)/api/$',
        views.fetch_document_relations, name='fetch-document-relations'),

    url(r'^(?P<relation_pk>\d+)/feedback-api/$',
        views.fetch_relation_feedback, name='fetch-relation-feedback'),

    url(r'^(?P<document_pk>\d+)/$',
        views.relation_task_home, name='task'),

)
