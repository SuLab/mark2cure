from django.contrib.auth.decorators import login_required
from django.shortcuts import get_object_or_404
from django.http import HttpResponse
from django.db import connection
from django.conf import settings

from ..document.models import Annotation, View
from ..task.models import Level, UserQuestRelationship

from .serializers import QuestSerializer, LeaderboardSerializer, NERGroupSerializer, TeamLeaderboardSerializer, DocumentRelationSerializer
from ..userprofile.models import Team
from ..common.models import Group
from ..analysis.models import Report
from ..task.models import Task
from ..task.entity_recognition.models import EntityRecognitionAnnotation
from ..task.relation.models import RelationAnnotation
from ..score.models import Point

from rest_framework.decorators import api_view
from rest_framework.response import Response

from itertools import chain
# from itertools import count
import networkx as nx
import datetime
import json


_attrs = dict(id='id', source='source', target='target', key='key')


def node_link_data(G, attrs=_attrs):
    multigraph = G.is_multigraph()
    id_ = attrs['id']

    source = attrs['source']
    target = attrs['target']

    # Allow 'key' to be omitted from attrs if the graph is not a multigraph.
    key = None if not multigraph else attrs['key']

    if len(set([source, target, key])) < 3:
        raise nx.NetworkXError('Attribute names are not unique.')

    data = {}
    data['directed'] = G.is_directed()
    data['multigraph'] = multigraph
    data['graph'] = list(G.graph.items())
    data['nodes'] = [dict(chain(G.node[n].items(), [(id_, n)])) for n in G]
    data['edges'] = [dict(chain(d.items(), [(source, u), (target, v), ('id', k)])) for u, v, k, d in G.edges_iter(keys=True, data=True)]  # N1, N2, IDX, ATTRS
    return data


def group_network(request, group_pk):
    group = get_object_or_404(Group, pk=group_pk)

    from ..analysis.tasks import generate_network
    G = generate_network(group.pk, spring_force=8)
    d = node_link_data(G)

    return HttpResponse(json.dumps(d), content_type='application/json')


@login_required
@api_view(['GET'])
def analysis_group_user(request, group_pk, user_pk=None):
    group = get_object_or_404(Group, pk=group_pk)

    response = []
    reports = group.report_set.filter(report_type=Report.AVERAGE).order_by('-created').all()
    user_id = int(user_pk) if user_pk else int(request.user.pk)

    for report in reports:
        df = report.dataframe
        df = df[df['user_id'] == user_id]

        if df.shape[0] > 0:
            row = df.iloc[0]
            response.append({
                'created': report.created,
                'f-score': row['f-score'],
                'pairings': row['pairings']})

    return Response(response)


@login_required
@api_view(['GET'])
def analysis_group(request, group_pk):
    group = get_object_or_404(Group, pk=group_pk)
    weighted = True

    response = []
    reports = group.report_set.filter(report_type=1).order_by('-created').all()
    for report in reports:
        df = report.dataframe

        if weighted:
            df['wf'] = df['pairings'] * df['f-score']
            response.append({
                'created': report.created,
                'f-score': df['wf'].sum() / df['pairings'].sum(),
                'pairings': df['pairings'].sum()})

        else:
            response.append({
                'created': report.created,
                'f-score': df['f-score'].mean(),
                'pairings': df['pairings'].sum()})

    return Response(response)


@api_view(['GET'])
def mark2cure_stats(request):
    return Response({
        'ner_annotations': EntityRecognitionAnnotation.objects.count(),
        're_annotations': RelationAnnotation.objects.count(),
    })


@login_required
@api_view(['GET'])
def user_task_stats(request):
    ner_level = Level.objects.filter(user=request.user, task_type='e').first()
    re_level = Level.objects.filter(user=request.user, task_type='r').first()

    return Response({
        'ner': ner_level.level if ner_level else 0,
        're': re_level.level if re_level else 0
    })


@login_required
@api_view(['GET'])
def ner_stats(request):
    return Response({
        'total_score': request.user.profile.score(task='entity_recognition'),
        'quests_completed': UserQuestRelationship.objects.filter(user=request.user, completed=True).count(),
        'papers_reviewed': View.objects.filter(user=request.user, completed=True, task_type='cr').count(),
        'annotations': Annotation.objects.filter(kind='e', view__user=request.user).count()
    })


@login_required
@api_view(['GET'])
def re_stats(request):
    return Response({
        'total_score': request.user.profile.score(task='relation'),
        'quests_completed': View.objects.filter(user=request.user, completed=True, task_type='ri').count(),
        'annotations': Annotation.objects.filter(kind='r', view__user=request.user).count()
    })


@login_required
@api_view(['GET'])
def ner_list_item_contributors(request, group_pk):
    group = get_object_or_404(Group, pk=group_pk)
    return Response([{'username': i[0], 'count': i[1]} for i in group.contributors()])


# @login_required
@api_view(['GET'])
def ner_list_item_quests(request, group_pk):
    group = get_object_or_404(Group, pk=group_pk)

    # we now allow users to see a group 'home page' for detailed information whether or
    # not they are logged in
    if request.user.is_authenticated():
        queryset = Task.objects.filter(kind=Task.QUEST, group=group).extra(select={
            "current_submissions_count": """
                SELECT COUNT(*) AS current_submissions_count
                FROM task_userquestrelationship
                WHERE (task_userquestrelationship.completed = 1
                    AND task_userquestrelationship.task_id = task_task.id)""",
            "user_completed": """
                SELECT COUNT(*) AS user_completed
                FROM task_userquestrelationship
                WHERE (task_userquestrelationship.completed = 1
                    AND task_userquestrelationship.user_id = %d
                    AND task_userquestrelationship.task_id = task_task.id)""" % (request.user.pk,)
        }).prefetch_related('documents')
    else:
        queryset = Task.objects.filter(kind=Task.QUEST, group=group).extra(select={
            "current_submissions_count": """
                SELECT COUNT(*) AS current_submissions_count
                FROM task_userquestrelationship
                WHERE (task_userquestrelationship.completed = 1
                    AND task_userquestrelationship.task_id = task_task.id)"""
        }).prefetch_related('documents')

    serializer = QuestSerializer(queryset, many=True, context={'user': request.user})
    return Response(serializer.data)


# @login_required
@api_view(['GET'])
def ner_list(request):
    queryset = Group.objects.exclude(stub='training').order_by('-order')
    serializer = NERGroupSerializer(queryset, many=True)
    return Response(serializer.data)


# @login_required
@api_view(['GET'])
def ner_list_item(request, group_pk):
    group = get_object_or_404(Group, pk=group_pk)

    # total annotation count here, plus anns


    # (TODO) this should not be hardcoded here; could not open file? Quick Fix. -JF
    # data should come from file in static/data "group_release_dates.txt"
    group_date_dict = {
        "CDG": {"invite": "2015.05.12", "public": "2015.05.21", "closed": "2015.07.29"},
        "alacrima": {"invite": "2015.05.21", "public": "2015.05.22", "closed": "2015.06.19"},
        "OGD": {"invite": "2015.05.29", "public": "2015.05.29", "closed": "2015.11.13"},
        "FBX": {"invite": "2015.06.25", "public": "2015.06.26", "closed": "2015.08.14"},
        "ost": {"invite": "2015.07.31", "public": "2015.08.07", "closed": "2016.04.03"},
        "mfold": {"invite": "2015.09.10", "public": "2015.09.11", "closed": "2015.11.19"},
        "eeyar": {"invite": "2015.11.04", "public": "2015.11.06", "closed": "2015.12.25"},
        "mitomis": {"invite": "2015.11.10", "public": "2015.11.11", "closed": "2016.03.04"},
        "ATGS": {"invite": "2015.12.28", "public": "2015.12.30", "closed": ""},
        "MATG": {"invite": "2016.02.24", "public": "2016.02.26", "closed": ""},
        "MATGS": {"invite": "2016.04.15", "public": "2016.04.15", "closed": ""},
        "training": {"invite": "2015.05.21", "public": "2015.05.21", "closed": ""},
        "HSPT1": {"invite": "2016.04.15", "public": "2016.04.15", "closed": ""},
    }
    try:
        start_date = group_date_dict[group.stub]['invite']
        end_date = group_date_dict[group.stub]['closed']
    except:
        start_date = ""
        end_date = ""

    return Response({
        "pk": 35,
        "name": group.name,
        "stub": group.stub,
        "description": group.description,
        "enabled": group.enabled,

        'document_count': group.document_count(),

        'total_contributors': 0,
        'percentage_complete': 0,
        "complete_percent": 70.66666666666667,
        'current_avg_f_score': 0,
        'start_date': start_date,
        'end_date': end_date
    })


@login_required
@api_view(['GET'])
def re_list(request):
    """ Returns the available relation tasks for a specific user
        Accessed through a JSON API endpoint
    """
    cmd_str = ""
    with open('mark2cure/api/commands/get-relations.sql', 'r') as f:
        cmd_str = f.read()

    # Start the DB Connection
    c = connection.cursor()

    c.execute('SET @user_work_max = {rel_work_size};'.format(rel_work_size=20))
    c.execute('SET @k_max = {completions};'.format(completions=settings.ENTITY_RECOGNITION_K))
    c.execute('SET @user_id = {user_id};'.format(user_id=request.user.pk))
    c.execute('SET @rel_ann_content_type_id = 56;')
    c.execute(cmd_str)

    queryset = [{'id': x[0],
                 'document_id': x[1],
                 'title': x[2],

                 'total_document_relationships': x[3],
                 'user_document_relationships': x[4],

                 'community_answered': x[5],
                 'community_completed': x[6],
                 'community_progress': x[7],

                 'user_completed': x[8],
                 'user_progress': x[9],
                 'user_answered': x[10],
                 'user_view_completed': x[11]} for x in c.fetchall()]

    # Close the connection
    c.close()

    serializer = DocumentRelationSerializer(queryset, many=True)
    return Response(serializer.data)


def users_with_score(days=30):
    today = datetime.datetime.now()
    since = today - datetime.timedelta(days=days)

    res = Point.objects.raw("""
        SELECT  ANY_VALUE(`score_point`.`id`) as `id`,
                SUM(score_point.amount) as score,
                `auth_user`.`username`,
                `auth_user`.`id`
        FROM `score_point`
        LEFT OUTER JOIN `auth_user`
            ON `auth_user`.`id` = `score_point`.`user_id`
        WHERE ( `score_point`.`created` > '{since}'
                AND `score_point`.`created` <= '{today}'
                AND `auth_user`.`id` NOT IN ({excluded_users}) )
        GROUP BY `auth_user`.`id` ORDER BY score DESC;""".format(
        since=since,
        today=today,
        excluded_users=', '.join('\'' + str(item) + '\'' for item in [5, 160]))
    )

    return [row for row in res if row.score is not None]


def get_annotated_teams(days=30):
    # (TODO) This could be smaller by only being UserProfiles that
    # we know are part of a Team
    users_queryset = users_with_score(days=days)

    teams = Team.objects.all()
    for team in teams:
        team_user_profile_pks = team.userprofile_set.values_list('pk', flat=True)
        team.score = sum(filter(None, [row.score for row in filter(lambda x: x.id in team_user_profile_pks, users_queryset)]))
    teams = list(teams)
    teams.sort(key=lambda x: x.score, reverse=True)
    return teams


@login_required
@api_view(['GET'])
def leaderboard_users(request, day_window):
    queryset = users_with_score(days=int(day_window))[:25]
    serializer = LeaderboardSerializer(queryset, many=True)
    return Response(serializer.data)


@login_required
@api_view(['GET'])
def leaderboard_teams(request, day_window):
    queryset = list(get_annotated_teams(days=int(day_window)))[:25]
    queryset = [team for team in queryset if team.score is not 0]
    serializer = TeamLeaderboardSerializer(queryset, many=True)
    return Response(serializer.data)

