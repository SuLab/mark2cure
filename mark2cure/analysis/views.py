from django.http import HttpResponse
from django.shortcuts import get_object_or_404

from django.contrib.auth.models import User
from mark2cure.common.models import Task, UserQuestRelationship

from django.contrib.auth.decorators import login_required
from mark2cure.common.bioc import *
from mark2cure.document.models import Document, Annotation

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
def document_bioc(request, doc_id):
    doc = get_object_or_404(Document, pk=doc_id)

    writer = BioCWriter()
    writer.collection = BioCCollection()
    collection = writer.collection
    collection.date = doc.updated.strftime("%Y%m%d")
    collection.source = doc.source

    d = BioCDocument()
    d.id = str(doc.document_id)

    passage_offset = 0
    for section in doc.available_sections():
        passage = BioCPassage()
        passage.put_infon('type', 'paragraph')
        passage.put_infon('section', section.get_kind_display())

        passage.offset = str(passage_offset)
        passage_offset += len(section.text)
        passage.text = section.text
        d.add_passage(passage)

        for ann in Annotation.objects.filter(view__section=section).all():
            annotation = BioCAnnotation()
            annotation.id = str(ann.pk)
            annotation.put_infon('type', str(ann.type))
            annotation.put_infon('user', str(ann.view.user.pk))

            location = BioCLocation()
            location.offset = str(passage_offset+ann.start)
            location.length = str(len(ann.text))
            annotation.add_location(location)

            annotation.text = ann.text

            passage.add_annotation(annotation)


    collection.add_document(d)

    return HttpResponse(writer, content_type='text/html')
