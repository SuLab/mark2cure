from django.conf.urls import patterns, url
from . import views

urlpatterns = patterns('',

    url(r'^$', views.home, name='home'),
    url(r'^beta/$', views.beta, name='beta'),
    url(r'^why-i-mark2cure/$', views.why_mark2cure, name='why-mark2cure'),

    url(r'^dashboard/$',
        views.dashboard, name='dashboard'),

    url(r'^support/$',
        views.support, name='support'),

    # Initial training for fresh signups
    url(r'^quest/(?P<quest_pk>\d+)/(?P<doc_idx>\d+)/$',
        views.quest_read_doc, name='quest-document'),

    url(r'^quest/(?P<quest_pk>\d+)/$',
        views.quest_read, name='quest-home'),

    url(r'^quest/(?P<quest_pk>\d+)/submit/$',
        views.quest_submit, name='quest-submit'),

    # REST Framework
    url(r'^quest/api/read/$',
        views.quest_list, name='quest-api'),
)
