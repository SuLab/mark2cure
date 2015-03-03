from django.conf.urls import patterns, url
from . import views


urlpatterns = patterns('',
    url(r'^users/training/$', views.users_training, name='users_training'),

    url(r'^(?P<doc_id>\d+).xml$',
        views.document_bioc, name='read-bioc'),
)
