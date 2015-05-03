from django.template.response import TemplateResponse
from django.shortcuts import get_object_or_404, redirect
from django.conf import settings
from django.views.decorators.http import require_http_methods

from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth.decorators import login_required
from django.contrib.messages import get_messages
from django.contrib import messages
from django.http import HttpResponse

from mark2cure.userprofile.models import UserProfile
from .models import Group, Task, UserQuestRelationship
from .serializers import QuestSerializer
from .forms import SupportMessageForm

from rest_framework.response import Response
from rest_framework.decorators import api_view

from brabeion import badges

import logging
import random
import os
logger = logging.getLogger(__name__)


@require_http_methods(['POST'])
def support(request):
    form = SupportMessageForm(data=request.POST)
    if form.is_valid():
        form.save()
        return HttpResponse(200)
    return HttpResponse(500)


def beta(request):
    return redirect('common:home')


def home(request):
    if request.user.is_authenticated():
        return redirect('common:dashboard')

    form = AuthenticationForm()
    quotes = ["To help others.", "In memory of my daughter who had Cystic Fibrosis.", "Rare disease dad!", "curiosity.", "This is needed.", "Goofing off productively.", "Community.", "Science!"]
    random.shuffle(quotes)
    return TemplateResponse(request, 'common/landing2.jade', {'form': form, 'quotes': quotes})


def why_mark2cure(request):
    query = UserProfile.objects.exclude(motivation='').order_by('?').values('motivation', 'user')
    return TemplateResponse(request, 'common/why-mark2cure.jade', {'profiles': query})


@login_required
def dashboard(request):
    if not request.user.profile.highest_level("skill").level == 7:
        return redirect('training:route')

    profile = request.user.profile
    tasks = Task.objects.filter(kind=Task.QUEST).all()
    for task in tasks:
        setattr(task, 'enabled', profile.highest_level('skill').level >= task.requires_qualification)
        setattr(task, 'completed', UserQuestRelationship.objects.filter(task=task, user=request.user, completed=True).exists())

    welcome = False
    storage = get_messages(request)
    for message in storage:
        if message.message == 'dashboard-unlock-success':
            welcome = True

    # Figure out state of the view for the user
    group = Group.objects.first()
    queryset = Task.objects.filter(kind=Task.QUEST, group=group).all()
    serializer = QuestSerializer(queryset, many=True, context={'user': request.user})

    user_completed = filter(lambda task: task['user']['completed'] is True, serializer.data)
    user_completed_count = len(user_completed)
    community_completed = filter(lambda task: task['progress']['completed'] is True, serializer.data)
    community_completed_count = len(community_completed)
    query_set_count = len(queryset)

    msg_footer = '<p class="text-center">Be sure to check your email and follow us on twitter (<a href="https://twitter.com/mark2cure">@Mark2Cure</a>) to be notified when we launch our next one.</p>'

    if user_completed_count == len(serializer.data):
        msg = '<p class="lead text-center">Congratulations! You have completed all quests available to you. Thank you for your participation in this experiment.</p>'

    elif user_completed_count >= 1 and community_completed_count == query_set_count - 1:
        msg = '<p class="lead text-center">Thank you very much for your participation in the first experiment. The Mark2Cure community has completed all the quests available.</p>'

    elif community_completed_count == query_set_count - 1:
        msg = '<p class="lead text-center">Thank you for joining Mark2Cure. The Mark2Cure community has completed all the quests available.</p>'

    else:
        msg = '<p class="lead text-center">Click on one of the quest numbers below to start the quest. Your contributions are important so complete as many quests as you can.</p>'

    messages.info(request, msg, extra_tags='safe alert-success')

    groups = []
    ctx = {
            'groups': groups,
            'tasks': tasks,
           'welcome': welcome,
           'profile': profile}
    return TemplateResponse(request, 'common/dashboard.jade', ctx)


@login_required
@api_view(['GET'])
def quest_list(request):
    group = Group.objects.first()
    queryset = Task.objects.filter(kind=Task.QUEST, group=group).all()
    serializer = QuestSerializer(queryset, many=True, context={'user': request.user})
    return Response(serializer.data)


def quest_prevent_duplicates(request, task):
    # Prevent a user from completing the same Quest multiple times
    if UserQuestRelationship.objects.filter(task=task, user=request.user, completed=True).exists():
        messages.warning(request, '<p class="lead text-center">We\'re sorry, but you can only do Quest {quest_pk} once.</p>'.format(quest_pk=task.name), extra_tags='safe alert-warning')
        return redirect('common:dashboard')


@login_required
def quest_read_doc(request, quest_pk, doc_idx):
    task = get_object_or_404(Task, pk=quest_pk)

    # Redirect if trying to access more documents
    # than this task contains
    if int(doc_idx) > len( task.remaining_documents() ):
        return redirect('common:quest-home', quest_pk=task.pk)

    # Confirm the user has started, but not completed this Task
    user_quest_relationship = task.user_relationship(request.user, False)
    if not user_quest_relationship:
        return redirect('common:quest-home', quest_pk=task.pk)

    task_doc_pks_completed = user_quest_relationship.completed_document_ids()
    if int(doc_idx) <= len( task_doc_pks_completed ):
        return redirect('common:quest-home', quest_pk=task.pk)

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
    random.shuffle( task_doc_uncompleted )

    document = task_doc_uncompleted[0]
    #task.create_views(document, request.user)

    ctx = {'task': task,
           'completed_doc_pks': task_doc_pks_completed,
           'uncompleted_docs': task_doc_uncompleted,
           'document': document}
    return TemplateResponse(request, 'common/quest.jade', ctx)


@login_required
def quest_read_doc_feedback(request, quest_pk, doc_idx):
    '''
        /quest/10/3/feedback/

        For now, we will not have a dedicated feedback page.
        This may be added in the future if desired, for now
        all feedback comes after ajax request to fetch the
        opponents BioC file.

        Logic for this decision: Would be 3 api requests
        User read, Opponent Read
                vs
        User read. User read, Opponent Read
    '''
    pass


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

        return redirect('common:dashboard')


@login_required
def quest_feedback(request, quest_pk):
    task = get_object_or_404(Task, pk=quest_pk)
    ctx = {'task':task}
    return TemplateResponse(request, 'common/quest-feedback.jade', ctx)


@login_required
def quest_read(request, quest_pk):
    task = get_object_or_404(Task, pk=quest_pk)

    # Check if user has pre-existing relationship with Quest
    user_quest_rel_queryset = UserQuestRelationship.objects.filter(task=task, user=request.user)
    #, completed=False)

    if user_quest_rel_queryset.exists():

        # User has viewed this Quest before and Completed it
        # so show them the feedback page
        if user_quest_rel_queryset.filter(completed=True).exists():
            return redirect('common:quest-feedback', quest_pk=task.pk)

        # Not using get_or_create b/c get occasionally returned multiple (unknown bug source)
        user_quest_relationship = user_quest_rel_queryset.first()
        task_doc_pks_completed = user_quest_relationship.completed_document_ids()

        # If there are no more documents to do, mark the Quest
        # as completed and go to dashboard
        task_doc_uncompleted = task.remaining_documents(task_doc_pks_completed)
        if len(task_doc_uncompleted) == 0:
            return quest_submit(request, task, True)

        next_doc_idx = len(task_doc_pks_completed)+1
        return redirect('common:quest-document', quest_pk=task.pk, doc_idx=next_doc_idx)

    else:
        # Create the User >> Quest relationship
        user_quest_rel = UserQuestRelationship.objects.create(task=task, user=request.user, completed=False)

        documents = list(task.documents.all())
        for document in documents:
            task.create_views(document, request.user)

        # Route the user to the right idx doc
        return redirect('common:quest-document', quest_pk=task.pk, doc_idx=1)

