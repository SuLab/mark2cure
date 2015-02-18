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
from .models import Task, UserQuestRelationship
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
    return TemplateResponse(request, 'common/index.jade', {'form': form, 'quotes': quotes})


def why_mark2cure(request):
    query = UserProfile.objects.exclude(motivation='').order_by('?').values('motivation', 'user')
    return TemplateResponse(request, 'common/why-mark2cure.jade', {'profiles': query})


@login_required
def dashboard(request):
    if request.user.profile.highest_level("skill").level <= 2:
        return redirect('training:index')

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
    queryset = Task.objects.filter(kind=Task.QUEST).all()
    serializer = QuestSerializer(queryset, many=True, context={'user': request.user})

    user_completed = filter(lambda task: task['user']['completed'] is True, serializer.data)
    user_completed_count = len(user_completed)
    community_completed = filter(lambda task: task['progress']['completed'] is True, serializer.data)
    community_completed_count = len(community_completed)
    query_set_count = len(queryset)
    msg_footer = '<p class="text-center">Be sure to check your email and follow us on twitter (<a href="https://twitter.com/mark2cure">@Mark2Cure</a>) to be notified when we launch our next one.</p>'

    if user_completed_count == len(serializer.data):
        msg = '<p class="lead text-center">Congratulations! You have completed all quests available to you. Thank you for your participation in this experiment.</p>'
        #messages.info(request, msg + msg_footer, extra_tags='safe alert-success')

    elif user_completed_count >= 1 and community_completed_count == query_set_count - 1:
        msg = '<p class="lead text-center">Thank you very much for your participation in the first experiment. The Mark2Cure community has completed all the quests available.</p>'
        #messages.info(request, msg + msg_footer, extra_tags='safe alert-info')

    elif community_completed_count == query_set_count - 1:
        msg = '<p class="lead text-center">Thank you for joining Mark2Cure. The Mark2Cure community has completed all the quests available.</p>'
        #messages.info(request, msg + msg_footer, extra_tags='safe alert-warning')

    else:
        msg = '<p class="lead text-center">Click on one of the quest numbers below to start the quest. Your contributions are important so complete as many quests as you can.</p>'
        #messages.info(request, msg, extra_tags='safe alert-success')

    msg = '<p class="lead text-center">The first Mark2Cure annotation campaign has just finished! Congratulations, nearly 9,000 annotation tasks were completed in 3 weeks thanks to our amazing community of Mark2Curators! While we analyze the collected data and prepare for the next annotation campaign, please do your part by helping us to find other Mark2Curators. Thank you! - Team Mark2Cure'
    messages.info(request, msg, extra_tags='safe alert-success')

    ctx = {'tasks': tasks,
           'welcome': welcome,
           'profile': profile}
    return TemplateResponse(request, 'common/dashboard.jade', ctx)


@login_required
@api_view(['GET'])
def quest_list(request):
    queryset = Task.objects.filter(kind=Task.QUEST, experiment=settings.EXPERIMENT).all()
    serializer = QuestSerializer(queryset, many=True, context={'user': request.user})
    return Response(serializer.data)


@login_required
def quest_read(request, quest_num):
    task = get_object_or_404(Task, pk=quest_num)
    #user_quest_rel, user_quest_rel_created = UserQuestRelationship.objects.get_or_create(task=task, user=request.user, completed=False)

    user_quest_rels = UserQuestRelationship.objects.filter(task=task, user=request.user, completed=False)
    user_quest_rel_created = False

    if user_quest_rels.exists():
        user_quest_rel = user_quest_rels.first()
    else:
        user_quest_rel = UserQuestRelationship.objects.create(task=task, user=request.user, completed=False)
        user_quest_rel_created = True

    task_doc_ids_completed = []

    if user_quest_rel_created:
        documents = list(task.documents.all())
        for document in documents:
            task.create_views(document, request.user)
    else:
        task_doc_ids_completed = list(set(user_quest_rel.views.filter(completed=True).values_list('section__document', flat=True)))
        documents = list(task.documents.exclude(pk__in=task_doc_ids_completed).all())

    random.shuffle(documents)

    if (request.method == 'POST' and user_quest_rel_created is False) or len(documents) == 0:
        # (TODO) Add validation check here at some point
        user_quest_rel.completed = True
        user_quest_rel.save()

        request.user.profile.rating.add(score=task.points, user=None, ip_address=os.urandom(7).encode('hex'))
        badges.possibly_award_badge("points_awarded", user=request.user)
        badges.possibly_award_badge("skill_awarded", user=request.user, level=task.provides_qualification)

        ctx = {'task': task}
        return TemplateResponse(request, 'common/quest-feedback.jade', ctx)

    user_quest_rel.save()
    ctx = {'task': task,
           'completed_docs': task_doc_ids_completed,
           'documents': documents}
    return TemplateResponse(request, 'common/quest.jade', ctx)

