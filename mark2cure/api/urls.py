from django.conf.urls import patterns, url
from . import views


urlpatterns = patterns('',
    # Analysis App

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
