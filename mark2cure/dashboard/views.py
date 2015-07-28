from django.contrib.auth.decorators import login_required
from django.template.response import TemplateResponse
from django.shortcuts import get_object_or_404
from django.contrib.auth.models import User
from django.http import HttpResponse

from ..common.models import Task, UserQuestRelationship, Group
from ..common.bioc import *

import csv


@login_required
def users_training(request):
    response = HttpResponse(content_type='text/html')
    writer = csv.writer(response)

    training = Task.objects.filter(kind=Task.TRAINING).all()
    writer.writerow(['user_pk', 'T1', 'T2', 'T3', 'E1_Q1', 'E2_Q2'])

    for user in User.objects.all():
        row = [user.pk]
        for t in training:
            row.append(UserQuestRelationship.objects.filter(task=t, user=user, completed=True).exists())

        row.append(UserQuestRelationship.objects.filter(task__pk=4, user=user, completed=True).exists())
        row.append(UserQuestRelationship.objects.filter(task__pk=11, user=user, completed=True).exists())
        writer.writerow(row)
    return response


@login_required
def group_list(request):
    ctx = {
        'groups': Group.objects.all()
    }
    return TemplateResponse(request, 'dashboard/group_list.jade', ctx)


@login_required
def group_read(request, pk):
    group = get_object_or_404(Group, pk=pk)
    ctx = {
        'group': group
    }
    return TemplateResponse(request, 'dashboard/group.jade', ctx)

