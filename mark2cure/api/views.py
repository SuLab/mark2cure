from django.shortcuts import render
from django.contrib.auth.decorators import login_required
from django.shortcuts import get_object_or_404
from django.http import HttpResponse, HttpResponseServerError

from rest_framework.response import Response
from rest_framework.decorators import api_view

from django.template.response import TemplateResponse

from mark2cure.api.serializers import QuestSerializer, UserProfileSerializer, GroupSerializer, TeamLeaderboardSerializer
from mark2cure.common.models import Group, Task, UserQuestRelationship
from mark2cure.common.formatter import bioc_writer, bioc_as_json
from mark2cure.userprofile.models import UserProfile, Team

from django.db.models import Sum
import datetime


@login_required
@api_view(['GET'])
def quest_group_list(request, group_pk):
    group = get_object_or_404(Group, pk=group_pk)
    queryset = Task.objects.filter(kind=Task.QUEST, group=group).all()
    serializer = QuestSerializer(queryset, many=True, context={'user': request.user})
    return Response(serializer.data)


def group_users_bioc(request, group_pk, format_type):
    group = get_object_or_404(Group, pk=group_pk)

    # When fetching via pubmed, include all user annotaitons
    writer = bioc_writer(request)

    for doc in group.get_documents():
        doc_bioc = doc.as_bioc_with_user_annotations()
        writer.collection.add_document(doc_bioc)

    if format_type == 'json':
        writer_json = bioc_as_json(writer)
        return HttpResponse(writer_json, content_type='application/json')
    else:
        return HttpResponse(writer, content_type='text/xml')


def group_pubtator_bioc(request, group_pk, format_type):
    group = get_object_or_404(Group, pk=group_pk)

    # When fetching via pubmed, include all user annotaitons
    writer = bioc_writer(request)

    for doc in group.get_documents()[:2]:
        doc_bioc = doc.as_bioc_with_user_annotations()
        writer.collection.add_document(doc_bioc)

    if format_type == 'json':
        writer_json = bioc_as_json(writer)
        return HttpResponse(writer_json, content_type='application/json')
    else:
        return HttpResponse(writer, content_type='text/xml')


@login_required
@api_view(['GET'])
def group_list(request):
    queryset = Group.objects.filter(enabled=True).order_by('-order').all()
    serializer = GroupSerializer(queryset, many=True)
    return Response(serializer.data)


def get_annotated_user_profiles(days=30, limit=25):
    '''
        content_type_id = 22 is a userprofile
    '''
    today = datetime.datetime.now()
    since = today - datetime.timedelta(days=days)

    return UserProfile.objects.exclude(pk__in=[5, 160]).extra(select = {
    "score" : """
        SELECT SUM(djangoratings_vote.score) AS score
        FROM djangoratings_vote
        WHERE (djangoratings_vote.content_type_id = 22
            AND djangoratings_vote.object_id = userprofile_userprofile.id
            AND djangoratings_vote.date_added > '%s'
            AND djangoratings_vote.date_added <= '%s')
        GROUP BY djangoratings_vote.object_id ORDER BY NULL""" % (since, today)
        }).prefetch_related('user').order_by("-score",)[:limit]


def get_annotated_teams(days=30, limit=25):
    '''
        content_type_id = 22 is a userprofile
    '''
    today = datetime.datetime.now()
    since = today - datetime.timedelta(days=days)

    # (TODO) This could be smaller by only being UserProfiles that
    # we know are part of a Team
    user_profiles = UserProfile.objects.exclude(pk__in=[5, 160]).extra(select = {
    "score" : """
        SELECT SUM(djangoratings_vote.score) AS score
        FROM djangoratings_vote
        WHERE (djangoratings_vote.content_type_id = 22
            AND djangoratings_vote.object_id = userprofile_userprofile.id
            AND djangoratings_vote.date_added > '%s'
            AND djangoratings_vote.date_added <= '%s')
        GROUP BY djangoratings_vote.object_id ORDER BY NULL""" % (since, today)
        }).order_by("-score",)

    teams = Team.objects.all()
    for team in teams:
        team_user_profile_pks = team.userprofile_set.values_list('pk', flat=True)
        team.score = sum(filter(None, user_profiles.filter(pk__in=team_user_profile_pks).values_list('score', flat=True)))
    teams = list(teams)
    teams.sort(key=lambda x: x.score, reverse=True)
    return teams


@login_required
@api_view(['GET'])
def leaderboard_users(request, day_window):
    queryset = get_annotated_user_profiles(days=int(day_window))
    serializer = UserProfileSerializer(queryset, many=True)
    return Response(serializer.data)


@login_required
@api_view(['GET'])
def leaderboard_teams(request, day_window):
    queryset = get_annotated_teams(days=int(day_window))
    serializer = TeamLeaderboardSerializer(queryset, many=True)
    return Response(serializer.data)

