from django.contrib.auth.decorators import login_required
from django.shortcuts import get_object_or_404
from django.http import HttpResponse

from .serializers import QuestSerializer, UserProfileSerializer, GroupSerializer, TeamLeaderboardSerializer
from ..common.formatter import bioc_writer, bioc_as_json
from ..userprofile.models import UserProfile, Team
from ..common.models import Group, Task

from rest_framework.decorators import api_view
from rest_framework.response import Response

import datetime


@login_required
@api_view(['GET'])
def quest_group_list(request, group_pk):
    group = get_object_or_404(Group, pk=group_pk)

    queryset = Task.objects.filter(kind=Task.QUEST, group=group).extra(select={
        "current_submissions_count": """
            SELECT COUNT(*) AS current_submissions_count
            FROM common_userquestrelationship
            WHERE (common_userquestrelationship.completed = 1
                AND common_userquestrelationship.task_id = common_task.id)""",
        "user_completed": """
            SELECT COUNT(*) AS user_completed
            FROM common_userquestrelationship
            WHERE (common_userquestrelationship.completed = 1
                AND common_userquestrelationship.user_id = %d
                AND common_userquestrelationship.task_id = common_task.id)""" % (request.user.pk,)
    }).prefetch_related('documents')

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

    for doc in group.get_documents():
        doc_bioc = doc.as_bioc_with_pubtator_annotations()
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


def userprofiles_with_score(days=30):
    today = datetime.datetime.now()
    since = today - datetime.timedelta(days=days)

    return UserProfile.objects.exclude(pk__in=[5, 160]).extra(select={
        "score": """
            SELECT SUM(djangoratings_vote.score) AS score
            FROM djangoratings_vote
            WHERE (djangoratings_vote.content_type_id = 22
                AND djangoratings_vote.object_id = userprofile_userprofile.id
                AND djangoratings_vote.date_added > '%s'
                AND djangoratings_vote.date_added <= '%s')
            GROUP BY djangoratings_vote.object_id ORDER BY NULL""" % (since, today)
    }).prefetch_related('user').order_by("-score",)


def get_annotated_teams(days=30):
    # (TODO) This could be smaller by only being UserProfiles that
    # we know are part of a Team
    userprofiles = userprofiles_with_score(days=days)

    teams = Team.objects.all()
    for team in teams:
        team_user_profile_pks = team.userprofile_set.values_list('pk', flat=True)
        team.score = sum(filter(None, userprofiles.filter(pk__in=team_user_profile_pks).values_list('score', flat=True)))
    teams = list(teams)
    teams.sort(key=lambda x: x.score, reverse=True)
    return teams


@login_required
@api_view(['GET'])
def leaderboard_users(request, day_window):
    queryset = list(userprofiles_with_score(days=int(day_window))[:25])
    queryset = [up for up in queryset if up.score is not None]
    serializer = UserProfileSerializer(queryset, many=True)
    return Response(serializer.data)


@login_required
@api_view(['GET'])
def leaderboard_teams(request, day_window):
    queryset = list(get_annotated_teams(days=int(day_window)))[:25]
    queryset = [team for team in queryset if team.score is not 0]
    serializer = TeamLeaderboardSerializer(queryset, many=True)
    return Response(serializer.data)

