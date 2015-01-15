from django.conf.urls import patterns, url
from . import views

urlpatterns = patterns('',

    url(r'^$', views.landing, name='landing'),
    url(r'^beta/$', views.home, name='home'),

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
