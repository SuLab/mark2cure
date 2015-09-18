from django.contrib.auth.decorators import login_required
from django.shortcuts import get_object_or_404
from django.http import HttpResponse

from ..common.bioc import BioCReader, BioCDocument, BioCPassage, BioCAnnotation, BioCLocation

from .serializers import QuestSerializer, UserProfileSerializer, GroupSerializer, TeamLeaderboardSerializer
from ..common.formatter import bioc_writer, bioc_as_json
from ..userprofile.models import UserProfile, Team
from ..common.models import Group, Task
from ..document.models import Section, Annotation

from rest_framework.decorators import api_view
from rest_framework.response import Response


import datetime
from networkx.readwrite import json_graph
import json


def network(request):

    from ..analysis.tasks import generate_network
    G = generate_network()

    d = json_graph.node_link_data(G)
    d['edges'] = d.pop('links')
    return HttpResponse(json.dumps(d), content_type='application/json')


@login_required
@api_view(['GET'])
def analysis_group_user(request, group_pk, user_pk=None):
    group = get_object_or_404(Group, pk=group_pk)

    if user_pk == None:
        user_pk = str(request.user.pk)

    response = []
    reports = group.report_set.filter(report_type=1).order_by('-created').all()
    for report in reports:
        df = report.dataframe
        df = df[df['user']==user_pk]
        if df.shape[0] > 0:
            row = df.iloc[0]
            response.append({
                'created': report.created,
                'f-score': row['f-score'],
                'pairings': row['pairings'] })

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
                'pairings': df['pairings'].sum() })

        else:
            response.append({
                'created': report.created,
                'f-score': df['f-score'].mean(),
                'pairings': df['pairings'].sum() })

    return Response(response)


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


#@api_view(['GET'])
def group_users_bioc(request, group_pk, format_type):
    '''
        Returns the BioC document for all user
        annotations accross the group
    '''
    content = False
    if not type(request) == dict:
        content = request.GET.get('content', False)

    # Fetch the group and all documents associated with the Group
    group = get_object_or_404(Group, pk=group_pk)
    document_pmids = group.get_documents().values_list('document_id', flat=True)

    # Fetch all the section pks and their text length
    sections = Section.objects.filter(
        document__document_id__in=document_pmids
    ).extra(select={'section_length': 'LENGTH(text)'
    }).values('pk', 'section_length')
    all_section_pks = [s['pk'] for s in sections]

    # Fetch all the actual annotations using
    annotations = Annotation.objects.filter(
            view__section__pk__in=all_section_pks
        ).values(
            'pk',
            'start',
            'text',
            'type',
            'view__user__pk',
            'view__section__pk',
            'view__section__document__document_id',
            ).all()
    annotations = list(annotations)

    # Provide the base of the response
    writer = bioc_writer(request)

    for doc_pmid in document_pmids:
        document_annotations = filter(lambda ann: ann['view__section__document__document_id'] == doc_pmid, annotations)

        document = BioCDocument()
        document.id = str(doc_pmid)

        passage_offset = 0

        section_pks = list(set([ann['view__section__pk'] for ann in document_annotations]))
        for section_pk in section_pks:
            # Add the Section to the Document
            passage = BioCPassage()
            passage.put_infon('id', str(section_pk))
            passage.offset = str(passage_offset)

            section_annotations = filter(lambda ann: ann['view__section__pk'] == section_pk, document_annotations)
            for ann in section_annotations:
                annotation = BioCAnnotation()

                annotation.id = str(ann.get('pk'))
                annotation.put_infon('user', str(ann.get('view__user__pk')))
                annotation.put_infon('type', ann.get('type'))

                location = BioCLocation()
                location.offset = str(passage_offset + ann.get('start'))
                location.length = str(len(ann.get('text')))
                annotation.add_location(location)
                annotation.text = ann.get('text')

                passage.add_annotation(annotation)

            section_results = filter(lambda section: section['pk'] == section_pk, sections)
            passage_offset += int(section_results[0]['section_length'])
            document.add_passage(passage)

        writer.collection.add_document(document)

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

