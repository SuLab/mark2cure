from django.conf.urls import url
from . import views


urlpatterns = [
    # Analysis App
    url(r'^network/(?P<group_pk>\d+)/$',
        views.group_network, name='group-network'),

    # Longitudinal user F in Group
    url(r'^analysis/group/(?P<group_pk>\d+)/user/(?P<user_pk>\d+)/$',
        views.analysis_group_user, name='analysis-group-user'),
    url(r'^analysis/group/(?P<group_pk>\d+)/user/$',
        views.analysis_group_user, name='analysis-group-user'),

    # Longitudinal Group F Avg
    url(r'^analysis/group/(?P<group_pk>\d+)/$',
        views.analysis_group, name='analysis-group'),

    # Tasks
    # - [Dashboard] Named Entity Recognition
    url(r'^ner/list/(?P<group_pk>\d+)/$',
        views.quest_group_list, name='quest-group-api'),

    url(r'^ner/list/$',
        views.group_list, name='groups-api'),

    # - [Dashboard] Relationship Extraction
    url(r're/list',
        views.relation_list, name='relations-api'),

    # - [Dashboard] User Scoreboard
    url(r'^leaderboard/users/(?P<day_window>\d+)/$',
        views.leaderboard_users, name='leaderboard-users'),

    url(r'^leaderboard/teams/(?P<day_window>\d+)/$',
        views.leaderboard_teams, name='leaderboard-teams')
]
