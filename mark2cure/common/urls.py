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
    url(r'^quest/(?P<quest_num>\w+)/$',
        views.quest_read, name='quest'),
    # REST Framework
    url(r'^quest/api/read/$',
        views.quest_list, name='quest-api'),
)
