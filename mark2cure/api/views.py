from django.shortcuts import render
from django.contrib.auth.decorators import login_required
from django.shortcuts import get_object_or_404

from rest_framework.response import Response
from rest_framework.decorators import api_view

from mark2cure.api.serializers import QuestSerializer, GroupSerializer
from mark2cure.common.models import Group, Task, UserQuestRelationship

# Create your views here.

@login_required
@api_view(['GET'])
def quest_group_list(request, group_pk):
    group = get_object_or_404(Group, pk=group_pk)
    queryset = Task.objects.filter(kind=Task.QUEST, group=group).all()
    serializer = QuestSerializer(queryset, many=True, context={'user': request.user})
    return Response(serializer.data)


@login_required
@api_view(['GET'])
def group_list(request):
    group = Group.objects.first()
    queryset = Group.objects.order_by('-pk').all()
    serializer = GroupSerializer(queryset, many=True)
    return Response(serializer.data)

