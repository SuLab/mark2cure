from django.conf.urls import patterns, url
from . import views


urlpatterns = patterns('',
    url(r'^users/training/$', views.users_training, name='users_training'),

    url(r'^bioc/all.xml$',
        views.document_bioc, name='read-bioc'),

    url(r'^bioc/(?P<doc_id>\d+).xml',
        views.document_bioc, name='read-bioc'),

    url(r'^bioc/(?P<doc_id>\d+).json',
        views.document_bioc_json, name='read-bioc-json'),

    url(r'^bioc/task/(?P<task_id>\d+)/(?P<doc_id>\d+)-opponent.xml',
        views.document_bioc_opponent, name='read-bioc-opponent'),
)

