from django.conf.urls import patterns, url
from . import views


urlpatterns = patterns('',
    # REST Framework
    url(r'^quest/(?P<group_pk>\d+)/$',
        views.quest_group_list, name='quest-group-api'),

    url(r'^groups/$',
        views.group_list, name='groups-api'),

)
