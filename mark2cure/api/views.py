from django.shortcuts import render
from django.contrib.auth.decorators import login_required
from django.shortcuts import get_object_or_404
from django.http import HttpResponse, HttpResponseServerError

from rest_framework.response import Response
from rest_framework.decorators import api_view

from django.template.response import TemplateResponse

from mark2cure.api.serializers import QuestSerializer, UserProfileSerializer, GroupSerializer
from mark2cure.common.models import Group, Task, UserQuestRelationship
from mark2cure.common.formatter import bioc_writer, bioc_as_json
from mark2cure.userprofile.models import UserProfile

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
    userprofiles = UserProfile.objects.all().extra(select = {
    "score" : """
        SELECT `djangoratings_vote`.`object_id`, SUM(`djangoratings_vote`.`score`) AS `score`
        FROM `djangoratings_vote` WHERE (`djangoratings_vote`.`content_type_id` = %s
        AND `djangoratings_vote`.`object_id` = %s
            AND `djangoratings_vote`.`date_added` > %s
            AND `djangoratings_vote`.`content_type_id` = %s
            AND `djangoratings_vote`.`date_added` <= %s)
            GROUP BY `djangoratings_vote`.`object_id` ORDER BY NULL""" % (user_id, 66, u'2015-05-11 20:46:42', 22, u'2015-06-10 20:46:42')
    }).order_by("-score",)
    '''


    # (TODO) Use Extra // http://timmyomahony.com/blog/filtering-annotations-django/
    sorted_up_dicts = []
    for up in UserProfile.objects.all():
        up_arr = up.rating.get_ratings().filter(content_type_id=22, date_added__lte=datetime.datetime.today(), date_added__gt=datetime.datetime.today()-datetime.timedelta(days=30)).values('object_id').annotate(score=Sum('score'))
        if len(up_arr) == 1:
            sorted_up_dicts.append(up_arr[0])

    sorted_up_dicts = sorted(sorted_up_dicts, key=lambda k: -k['score'])
    return sorted_up_dicts[:limit]


@login_required
@api_view(['GET'])
def leaderboard_users(request):
    #userprofiles = get_annotated_user_profiles(30)
    queryset = UserProfile.objects.exclude(user__username__in=['CheckingOnTesters', 'roarwithphoebe']).order_by('-rating_score').all()[:25]
    serializer = UserProfileSerializer(queryset, many=True)
    return Response(serializer.data)
