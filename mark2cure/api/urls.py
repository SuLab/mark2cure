from django.conf.urls import patterns, url
from . import views


urlpatterns = patterns('',
    # REST Framework
    url(r'^quest/(?P<group_pk>\d+)/$',
        views.quest_group_list, name='quest-group-api'),

    # BioC File for M2C annotations
    url(r'^group/(?P<group_pk>\d+)/user_annotations.(?P<format_type>\w+)$',
        views.quest_users_bioc, name='quest-users-bioc'),

    url(r'^groups/$',
        views.group_list, name='groups-api'),

    url(r'^leaderboard/users/$',
        views.leaderboard_users, name='leaderboard-users'),

)
