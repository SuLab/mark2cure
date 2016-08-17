from django.contrib.auth.decorators import login_required
from django.shortcuts import get_object_or_404
from django.http import HttpResponse

from .serializers import QuestSerializer, LeaderboardSerializer, GroupSerializer, TeamLeaderboardSerializer, DocumentRelationSerializer
from ..common.formatter import bioc_writer, bioc_as_json, clean_df, apply_annotations
from ..userprofile.models import Team
from ..common.models import Document, Group
from ..task.models import Task
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

    if user_pk is None:
        user_pk = str(request.user.pk)

    response = []
    reports = group.report_set.filter(report_type=1).order_by('-created').all()
    for report in reports:
        df = report.dataframe
        df = df[df['user'] == user_pk]
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


# @login_required
@api_view(['GET'])
def quest_group_list(request, group_pk):
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


@login_required
@api_view(['GET'])
def relation_list(request):

    queryset = Document.objects.filter(relationgroup__stub='alacrima').extra(select={
        "current_completed_count": """
            SELECT COUNT(*) AS user_completed_count
            FROM document_view

            INNER JOIN `document_section` ON ( document_view.section_id = document_section.id)

            WHERE (document_view.task_type = 'ri'
                AND document_view.completed = 1
                AND document_section.document_id = document_document.id)""",

        "task_count": """
            SELECT COUNT(*) AS task_count
            FROM relation_relation
            WHERE relation_relation.document_id = document_document.id""",

        "concepts": """
            SELECT COUNT(*) AS task_count
            FROM relation_relation
            WHERE relation_relation.document_id = document_document.id""",

        # How many times has this Document Task View been completed (should only be 0 or 1)
        "user_completed_count": """
            SELECT COUNT(*) AS user_completed_count
            FROM document_view

            INNER JOIN `document_section` ON ( document_view.section_id = document_section.id)

            WHERE (document_view.task_type = 'ri'
                AND document_view.completed = 1
                AND document_view.user_id = %d
        AND document_section.document_id = document_document.id)""" % (request.user.pk,)
    })

    queryset = [x for x in queryset if x.task_count > 0]
    queryset = [x for x in queryset if x.task_count < 20]
    queryset = [x for x in queryset if x.user_completed_count == 0]

    serializer = DocumentRelationSerializer(queryset, many=True, context={'user': request.user})
    return Response(serializer.data)


def group_bioc(request, group_pk, selection_type, format_type):
    '''
        Returns the BioC document for all annotations accross the group
    '''

    # Fetch the group and all documents associated with the Group
    group = get_object_or_404(Group, pk=group_pk)

    writer = bioc_writer(request)

    for doc in group.get_documents():
        doc_writer = doc.as_writer()
        if selection_type == 'user':
            doc_df = doc.as_er_df_with_user_annotations()
        else:
            doc_df = doc.as_er_df_with_pubtator_annotations()

        # Protection isn't needed b/c this is the raw output for analysis.
        doc_df = clean_df(doc_df, overlap_protection=False)

        # convert DF table into BioC Document
        doc_writer = apply_annotations(doc_writer, doc_df)

        writer.collection.add_document(doc_writer.collection.documents[0])

    if format_type == 'json':
        writer_json = bioc_as_json(writer)
        return HttpResponse(writer_json, content_type='application/json')
    else:
        return HttpResponse(writer, content_type='text/xml')


# @login_required
@api_view(['GET'])
def group_list(request):
    queryset = Group.objects.exclude(stub='training').order_by('-order')
    serializer = GroupSerializer(queryset, many=True)
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

