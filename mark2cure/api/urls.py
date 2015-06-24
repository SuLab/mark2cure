from django.conf.urls import patterns, url
from . import views


urlpatterns = patterns('',
    # REST Framework
    url(r'^quest/(?P<group_pk>\d+)/$',
        views.quest_group_list, name='quest-group-api'),

    # BioC File for M2C annotations by Group
    url(r'^group/(?P<group_pk>\d+)/user_annotations.(?P<format_type>\w+)$',
        views.group_users_bioc, name='group-users-bioc'),

    # BioC File for Pubtator annotations by Group
    url(r'^group/(?P<group_pk>\d+)/pubtator_annotations.(?P<format_type>\w+)$',
        views.group_pubtator_bioc, name='group-pubtator-bioc'),

    url(r'^groups/$',
        views.group_list, name='groups-api'),

    url(r'^leaderboard/users/(?P<day_window>\d+)/$',
        views.leaderboard_users, name='leaderboard-users'),

    url(r'^leaderboard/teams/(?P<day_window>\d+)/$',
        views.leaderboard_teams, name='leaderboard-teams'),

)
