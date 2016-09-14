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

    # REST Framework
    url(r'^quest/(?P<group_pk>\d+)/$',
        views.quest_group_list, name='quest-group-api'),

    url(r'^groups/$',
        views.group_list, name='groups-api'),

    url(r'relationships',
        views.relation_list, name='relations-api'),

    url(r'^leaderboard/users/(?P<day_window>\d+)/$',
        views.leaderboard_users, name='leaderboard-users'),

    url(r'^leaderboard/teams/(?P<day_window>\d+)/$',
        views.leaderboard_teams, name='leaderboard-teams')
]
