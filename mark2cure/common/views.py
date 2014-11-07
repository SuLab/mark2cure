from django.template import RequestContext
from django.shortcuts import get_object_or_404, render_to_response, redirect
from django.views.decorators.http import require_http_methods

from django.contrib.auth.models import User
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth.decorators import login_required
from django.contrib.auth import authenticate, login, logout
from django.contrib.messages import get_messages
from django.contrib import messages

from django.http import HttpResponse
from django.conf import settings
from django.core.mail import send_mail
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger

from mark2cure.document.models import Document, View, Annotation
from mark2cure.common.models import Task, UserQuestRelationship
from mark2cure.common.serializers import QuestSerializer
from rest_framework import viewsets, generics
from rest_framework.response import Response
from rest_framework.decorators import api_view

from brabeion import badges
from brabeion.models import BadgeAward
#from zinnia.models import Entry

from datetime import datetime, timedelta
import math
import random
import os
import logging
logger = logging.getLogger(__name__)


def signup_home(request):
    return render_to_response('common/landing.jade', {}, context_instance=RequestContext(request))


def home(request):

    if request.user.is_authenticated():
      return redirect('mark2cure.common.views.dashboard')

    form = AuthenticationForm()
    return render_to_response('common/index.jade',
                              {'form': form},
                              context_instance=RequestContext(request))


def introduction(request):
    return render_to_response('training/basics.jade',
                              {}, context_instance=RequestContext(request))


def training_read(request):
    tasks = []
    if request.user.is_authenticated():
        profile = request.user.profile

    tasks = Task.objects.filter(kind=Task.TRAINING).all()
    for task in tasks:
        setattr(task, 'enabled', profile.highest_level('skill').level >= task.requires_qualification if request.user.is_authenticated() else False)
        setattr(task, 'completed', UserQuestRelationship.objects.filter(task=task, user=request.user, completed=True).exists() if request.user.is_authenticated() else False)

        if task.pk == 1:
            setattr(task, 'enabled', True)

    return render_to_response('training/training_read.jade',
                              {'tasks': tasks},
                              context_instance=RequestContext(request))


def training_one(request, step_num):
    if step_num == 'complete':
        return render_to_response(
            'training/intro-1/complete.jade',
            {}, context_instance=RequestContext(request))

    step_num = int(step_num)
    next_ = step_num+1
    if step_num == 1:
        header1 = "Let's start by marking diseases"
        header2 = "Mark all the disease terms in the sentence below."
        paragraph = "Does choice of insulin regimen really matter in the management of diabetes ?"
        answers = [{'text': 'diabetes', 'start': 66}]

    if step_num == 2:
        header1 = "Sometimes you will see multiple instances of the same disease - Be sure to mark them all!"
        header2 = "Mark all the disease terms in the sentence below."
        paragraph = "To assess the management of diabetes , we reviewed records of 20 diabetes patients ."
        answers = [{'text': 'diabetes', 'start': 28}, {'text': 'diabetes', 'start': 65}]

    if step_num == 3:
        header1 = "Sometimes the disease is described by a conjuction of several words. Mark these disease conjunctions as a single span of text."
        header2 = "Mark all the disease terms in the sentence below."
        paragraph = "Of the 20 diabetes patients , 17 had type 2 diabetes mellitus ."
        answers = [{'text': 'diabetes', 'start': 10}, {'text': 'type 2 diabetes mellitus', 'start': 37}];

    if step_num == 4:
        header1 = "Sometimes the disease conjunctions are separated by words like 'and/or'. Decide if there are two distinct diseases to highlight, or if it is a single disease conjunction."
        header2 = "Mark all the disease terms in the sentence below."
        paragraph = "The remaining 3 had inherited and/or type I diabetes mellitus ."
        answers = [{'text': 'inherited and/or type I diabetes mellitus', 'start': 20}];

    if step_num == 5:
        header1 = "Sometimes different diseases are discussed. Mark all the diseases below."
        header2 = "Remember to mark disease conjunctions as spans and to mark different diseases separately."
        paragraph = "Of the 20 patients , 10 patients were also diagnosed with heart disease or rheumatoid arthritis ."
        answers = [{'text': 'heart disease', 'start': 58}, {'text': 'rheumatoid arthritis', 'start': 75}];

    if step_num == 6:
        header1 = "Sometimes the disease terms are abbreviated. Mark all instances of disease abbreviations."
        header2 = "Mark all disease terms in the sentence below."
        paragraph = "We will discuss the effect of different insulin regimen on type 2 diabetes mellitus patients (ie- T2DM patients) ..."
        answers = [{'text': 'type 2 diabetes mellitus', 'start': 59}, {'text': 'T2DM', 'start': 98}];

    if step_num == 7:
        header1 = "Practice what you've learned so far..."
        header2 = "Try marking the disease and disease abbreviations in this phrase now!"
        paragraph = "... with or without rheumatoid arthritis ( RA ) or heart disease ( HD ) ."
        answers = [{'text':'rheumatoid arthritis', 'start': 20}, {'text': 'RA', 'start': 43}, {'text': 'heart disease', 'start': 51}, {'text': 'HD', 'start': 67}];

    if step_num == 8:
        header1 = "Now let's mark some symptoms."
        header2 = "Symptoms are the physical manifestations of the disease. Mark the symptoms in the sentence below."
        paragraph = "In particular, we will examine the effects of these insulin regimen on the symptoms of the diseases. We will focus on some common symptoms such as fatigue as well as ..."
        answers = [{'text': 'fatigue', 'start': 147}];

    if step_num == 9:
        header1 = "Sometimes a single symptom is described using more than one word."
        header2 = "Mark these symptoms as single spans of text. Mark the symptom in the text below."
        paragraph = "... as well as frequent urination . Another crucial symptom ..."
        answers = [{'text': 'frequent urination', 'start': 15}];

    if step_num == 10:
        header1 = "Sometimes a single symptom is described with a long block of text and may have joining terms such as 'and' or 'or'."
        header2 = "In that case, highlight the entire symptom as a single span of text. Finish Training #1 by marking the symptom in the text below."
        paragraph = "Another crucial symptom that we will investigate includes tingling and/or numbness in the hands or feet ."
        answers = [{'text': 'tingling and/or numbness in the hands or feet', 'start': 58}];
        next_ = 'complete'

    return render_to_response(
        'training/intro-1/read.jade',
        {'training_num': 1,
         'step_num': step_num,
         'header1': header1,
         'header2': header2,
         'paragraph': paragraph,
         'answers': answers,
         'next': next_},
        context_instance=RequestContext(request))

@login_required
def training_two(request, step_num):
    if step_num == 'complete':
        if request.user.is_authenticated():
            task = Task.objects.get(pk=2)
            UserQuestRelationship.objects.create(task=task, user=request.user, completed=True)
            request.user.profile.rating.add(score=task.points, user=None, ip_address=os.urandom(7).encode('hex'))
            badges.possibly_award_badge("points_awarded", user=request.user)


        badges.possibly_award_badge("skill_awarded", user=request.user, level=3)
        return redirect('mark2cure.common.views.training_read')

    return render_to_response(
        'training/intro-2/step-{step_num}.jade'.format(step_num=step_num),
        {'step_num': step_num},
        context_instance=RequestContext(request))


@login_required
def training_three(request):
    if request.method == 'POST':
        task = Task.objects.get(pk=3)
        UserQuestRelationship.objects.create(task=task, user=request.user, completed=True)
        request.user.profile.rating.add(score=task.points, user=None, ip_address=os.urandom(7).encode('hex'))
        badges.possibly_award_badge("points_awarded", user=request.user)
        badges.possibly_award_badge("skill_awarded", user=request.user, level=4)
        messages.success(request, 'dashboard-unlock-success')

        return redirect('mark2cure.common.views.dashboard')

    else:
        return render_to_response(
            'training/intro-3/read.jade', {},
            context_instance=RequestContext(request))


@login_required
def dashboard(request):
    if request.user.profile.highest_level("skill").level <= 2:
       return redirect('mark2cure.common.views.training_read')

    profile = request.user.profile
    # posts = Entry.objects.all()[:3]
    posts = []
    tasks = Task.objects.filter(kind=Task.QUEST).all()
    for task in tasks:
        setattr(task, 'enabled', profile.highest_level('skill').level >= task.requires_qualification)
        setattr(task, 'completed', UserQuestRelationship.objects.filter(task=task, user=request.user, completed=True).exists())

    welcome = False
    storage = get_messages(request)
    for message in storage:
        if message.message == 'dashboard-unlock-success':
            welcome = True

    return render_to_response('common/dashboard.jade',
                              {'posts': posts,
                               'tasks': tasks,
                               'welcome': welcome,
                               'profile': profile},
                              context_instance=RequestContext(request))


@api_view(['GET'])
def quest_list(request):
    queryset = Task.objects.filter(kind=Task.QUEST).all()
    serializer = QuestSerializer(queryset, many=True, context={'user': request.user})
    return Response(serializer.data)


@login_required
def quest_read(request, quest_num):
    task = get_object_or_404(Task, pk=quest_num)
    user_quest_rel, user_quest_rel_created = UserQuestRelationship.objects.get_or_create(task=task, user=request.user, completed=False)

    if user_quest_rel_created:
        documents = list(task.documents.all())
        for document in documents:
            task.create_views(document, request.user)
    else:
        task_doc_ids_completed = list(set(user_quest_rel.views.filter(completed=True).values_list('section__document', flat=True)))
        documents = list(task.documents.exclude(pk__in=task_doc_ids_completed).all())

    random.shuffle(documents)

    if request.method == 'POST' and user_quest_rel_created is False:
        # (TODO) Add validation check here at some point
        user_quest_rel.completed = True
        user_quest_rel.save()

        request.user.profile.rating.add(score=task.points, user=None, ip_address=os.urandom(7).encode('hex'))
        badges.possibly_award_badge("points_awarded", user=request.user)
        badges.possibly_award_badge("skill_awarded", user=request.user, level=task.provides_qualification)
        return redirect('mark2cure.common.views.dashboard')

    user_quest_rel.save()
    return render_to_response('common/quest.jade',
                              {'task': task,
                               'documents': documents},
                              context_instance=RequestContext(request))


@login_required
def profile_survey(request):
    user_profile = request.user.userprofile

    if request.method == 'POST':
        form = ProfileSurveyForm(request.POST, instance = user_profile)
        if form.is_valid():
            form.save()
            return redirect('mark2cure.common.views.library')

    else:
        form = ProfileSurveyForm(instance = user_profile)

    return render_to_response('common/survey.jade',
                              {'user_profile': user_profile, 'form': form },
                              context_instance=RequestContext(request))


@require_http_methods(["POST"])
@login_required
def message(request):
    form = MessageForm(request.POST)
    if form.is_valid():
      message = form.save(commit=False)
      message.user = request.user
      message.save()
      return HttpResponse("Success")

    return HttpResponse('Unauthorized', status=401)


@require_http_methods(["POST"])
@login_required
def survey(request):
    for k, v in request.POST.iteritems():
        # print k, v
        if(k != "csrfmiddlewaretoken"):
          sf = SurveyFeedback(question = k, response = v, user = request.user)
          sf.save()

    return HttpResponse("Success")



