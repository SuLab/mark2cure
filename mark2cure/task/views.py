from django.http import HttpResponse, HttpResponseServerError
from django.views.decorators.http import require_http_methods
from django.contrib.auth.decorators import login_required
from django.shortcuts import get_object_or_404, redirect
from django.template.response import TemplateResponse
from django.views.decorators.cache import cache_page
from django.contrib.auth.models import User
from django.contrib import messages

from ..common.formatter import bioc_as_json, apply_bioc_annotations
from ..document.models import Document
from .models import Task, UserQuestRelationship

from brabeion import badges

import logging
import random
import os
logger = logging.getLogger(__name__)


@login_required
def quest_read_doc(request, quest_pk, doc_idx):
    task = get_object_or_404(Task, pk=quest_pk)

    # Redirect if trying to access more documents
    # than this task contains
    if int(doc_idx) > task.remaining_documents_count():
        return redirect('task:quest-home', quest_pk=task.pk)

    # Confirm the user has started, but not completed this Task
    user_quest_relationship = task.user_relationship(request.user, False)
    if not user_quest_relationship:
        return redirect('task:quest-home', quest_pk=task.pk)

    task_doc_pks_completed = user_quest_relationship.completed_document_ids()
    if int(doc_idx) <= len(task_doc_pks_completed):
        return redirect('task:quest-home', quest_pk=task.pk)

    '''
    # only take the one that has views
    user_quest_rel = task.userquestrelationship_set.filter(user=user, completed=False).latest()
    user_quest_rel_views = user_quest_rel.views

    # If a completed view doesn't exist, redirect them to complete the document submission
    if not user_quest_rel_views.filter(section__document=doc, completed=True).exists():
        return redirect('document:read', task.pk, doc.pk)
    '''

    # Fetch available documents
    task_doc_uncompleted = task.remaining_documents(task_doc_pks_completed)
    random.shuffle(task_doc_uncompleted)

    document = task_doc_uncompleted[0]
    # task.create_views(document, request.user)

    ctx = {'task': task,
           'completed_doc_pks': task_doc_pks_completed,
           'uncompleted_docs': task_doc_uncompleted,
           'document': document}
    return TemplateResponse(request, 'entity_recognition/quest.jade', ctx)


@login_required
def quest_read_doc_results_bioc(request, quest_pk, doc_pk, user_pk, format_type):
    # (TODO) remove task, why was this passed?
    task = get_object_or_404(Task, pk=quest_pk)
    document = get_object_or_404(Document, pk=doc_pk)
    user = get_object_or_404(User, pk=user_pk)

    # BioC Writer Response that will serve all partner comparison information
    writer = document.as_writer()
    writer = apply_bioc_annotations(writer, user)

    writer.collection.put_infon('partner', str(user.username))
    writer.collection.put_infon('partner_level', str(user.userprofile.highest_level().name))

    if format_type == 'json':
        writer_json = bioc_as_json(writer)
        return HttpResponse(writer_json, content_type='application/json')
    else:
        return HttpResponse(writer, content_type='text/xml')


@login_required
def quest_read_doc_results(request, quest_pk, doc_idx):
    '''
        Allows player to revist result page to look at comparision
    '''
    task = get_object_or_404(Task, pk=quest_pk)
    # Get the UQR. Completed=False b/c they may want to view results
    # of previous doc_idx's before finishing the rest of the document
    user_quest_relationship = task.user_relationship(request.user, True)

    # The completed document at this index
    # Using Abstracts but this [1,1,2,2,3,3,4,4,5,5] assumption
    # is extremely fragile
    relevant_views = user_quest_relationship.views.filter(section__kind='a')

    if relevant_views.filter(opponent__isnull=False).exists():
        relevant_view = relevant_views.filter(opponent__isnull=False)[int(doc_idx)]
        opponent = relevant_view.opponent.user
    else:
        # Novel annotations for the player (unpaired)
        relevant_view = relevant_views.all()[int(doc_idx)]
        opponent = None

    ctx = {'task': task,
           'opponent': opponent,
           'document': relevant_view.section.document}
    return TemplateResponse(request, 'entity_recognition/quest-results.jade', ctx)



@login_required
@require_http_methods(['POST'])
def document_quest_submit(request, quest_pk, document_pk):
    task = get_object_or_404(Task, pk=quest_pk)
    user_quest_relationship = task.user_relationship(request.user, False)

    if not user_quest_relationship:
        return HttpResponseServerError()

    for view in user_quest_relationship.views.filter(section__document__pk=document_pk):
        view.completed = True
        view.save()

    return HttpResponse(200)


@login_required
def quest_submit(request, task, bypass_post=False):
    # (TODO) Add validation check here at some point

    if request.POST or bypass_post:
        user_quest_relationship = task.user_relationship(request.user, False)

        if not user_quest_relationship.completed:
            request.user.profile.rating.add(score=task.points, user=None, ip_address=os.urandom(7).encode('hex'))
            badges.possibly_award_badge("points_awarded", user=request.user)
            badges.possibly_award_badge("skill_awarded", user=request.user, level=task.provides_qualification)

        user_quest_relationship.completed = True
        user_quest_relationship.save()


@login_required
def quest_feedback(request, quest_pk):
    task = get_object_or_404(Task, pk=quest_pk)
    ctx = {'task': task}
    return TemplateResponse(request, 'entity_recognition/quest-feedback.jade', ctx)


@login_required
def quest_read(request, quest_pk):
    task = get_object_or_404(Task, pk=quest_pk)

    # Check if user has pre-existing relationship with Quest
    user_quest_rel_queryset = UserQuestRelationship.objects.filter(task=task, user=request.user)

    if user_quest_rel_queryset.exists():

        # User has viewed this Quest before and Completed it
        # so show them the feedback page
        if user_quest_rel_queryset.filter(completed=True).exists():
            return redirect('task:quest-feedback', quest_pk=task.pk)

        # Not using get_or_create b/c get occasionally returned multiple (unknown bug source)
        user_quest_relationship = user_quest_rel_queryset.first()
        task_doc_pks_completed = user_quest_relationship.completed_document_ids()

        # If there are no more documents to do, mark the Quest
        # as completed and go to dashboard
        task_doc_uncompleted = task.remaining_documents(task_doc_pks_completed)
        if len(task_doc_uncompleted) == 0:
            quest_submit(request, task, True)
            return redirect('task:quest-feedback', quest_pk=task.pk)

        next_doc_idx = len(task_doc_pks_completed) + 1
        return redirect('task:quest-document', quest_pk=task.pk, doc_idx=next_doc_idx)

    else:
        # Create the User >> Quest relationship
        UserQuestRelationship.objects.create(task=task, user=request.user, completed=False)

        documents = list(task.documents.all())
        for document in documents:
            task.create_views(document, request.user)

        # Route the user to the right idx doc
        return redirect('task:quest-document', quest_pk=task.pk, doc_idx=1)

