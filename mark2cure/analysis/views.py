from django.http import HttpResponse
from django.shortcuts import get_object_or_404

from django.contrib.auth.models import User
from mark2cure.common.models import Task, UserQuestRelationship

from django.contrib.auth.decorators import login_required
from mark2cure.common.bioc import *
from mark2cure.document.models import Document, Annotation
from mark2cure.analysis.utils import apply_bioc_documents
from mark2cure.document.utils import select_best_opponent

from mark2cure.common.models import Task

import csv
import xmltodict
import json


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
def document_bioc(request, doc_id):
    query_set = Document.objects.filter(pk=doc_id)

    writer = BioCWriter()
    writer.collection = BioCCollection()
    collection = writer.collection
    collection.date = query_set.first().updated.strftime("%Y%m%d")
    collection.source = query_set.first().source

    apply_bioc_documents(query_set.all(), collection)

    return HttpResponse(writer, content_type='text/xml')


@login_required
def document_bioc_json(request, doc_id):
    query_set = Document.objects.filter(pk=doc_id)

    writer = BioCWriter()
    writer.collection = BioCCollection()
    collection = writer.collection
    collection.date = query_set.first().updated.strftime("%Y%m%d")
    collection.source = query_set.first().source

    apply_bioc_documents(query_set.all(), collection)

    o = xmltodict.parse(writer.__str__())
    print json.dumps(o)
    #json.dumps(o) # '{"e": {"a": ["text", "text"]}}''}'
    return HttpResponse(json.dumps(o), content_type='application/json')


@login_required
def document_bioc_opponent(request, task_id, doc_id):
    task = get_object_or_404(Task, pk=task_id)
    # (TODO) follow common_document_quest_relationship
    query_set = Document.objects.filter(pk=doc_id)

    opponent = select_best_opponent(task, query_set.first(), request.user)

    writer = BioCWriter()
    writer.collection = BioCCollection()
    collection = writer.collection
    collection.date = query_set.first().updated.strftime("%Y%m%d")
    collection.source = query_set.first().source

    apply_bioc_documents(query_set.all(), collection, opponent)

    return HttpResponse(writer, content_type='text/xml')
