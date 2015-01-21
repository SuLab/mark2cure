from django.shortcuts import render
from django.http import HttpResponse

from django.contrib.auth.models import User
from mark2cure.common.models import Task, UserQuestRelationship

import csv

def users_training(request):
    response = HttpResponse(content_type='text/html')
    #response['Content-Disposition'] = 'attachment; filename="users_training_completion.csv"'
    writer = csv.writer(response)

    training = Task.objects.filter(kind=Task.TRAINING).all()
    writer.writerow(['username', 'email', 'T1', 'T2', 'T3', 'E1_Q1', 'E2_Q2'])

    for user in User.objects.all():
        row = [user.username, user.email]
        for t in training:
            row.append(UserQuestRelationship.objects.filter(task=t, user=user, completed=True).exists())

        row.append(UserQuestRelationship.objects.filter(task__pk=4, user=user, completed=True).exists())
        row.append(UserQuestRelationship.objects.filter(task__pk=11, user=user, completed=True).exists())
        writer.writerow(row)
    return response
