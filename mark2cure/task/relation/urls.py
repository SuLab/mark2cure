from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^(?P<document_pk>\d+)/analysis/(?P<relation_pk>\d+)/$',
        views.document_analysis, name='document-analysis-specific'),
    url(r'^(?P<document_pk>\d+)/analysis/$',
        views.document_analysis, name='document-analysis'),


    url(r'^(?P<document_pk>\d+)/api/$',
        views.re_task_relationships_list, name='fetch-document-relations'),

    url(r'^(?P<document_pk>\d+)/(?P<relation_pk>\d+)/submit/$',
        views.re_task_relationship_submit, name='re-task-relationship-submit'),

    url(r'^(?P<document_pk>\d+)/submit/$',
        views.re_task_submit, name='re-task-submit'),

    url(r'^(?P<document_pk>\d+)/$',
        views.re_task, name='re-task'),

]
