from django.shortcuts import render
from django.template.response import TemplateResponse
from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse
from django.shortcuts import redirect
from django.contrib import messages

from mark2cure.common.models import Task, UserQuestRelationship

from brabeion import badges
from brabeion.models import BadgeAward

import os

def introduction(request):
    return TemplateResponse(request, 'training/basics.jade')


def read(request):
    tasks = []
    if request.user.is_authenticated():
        profile = request.user.profile

    tasks = Task.objects.filter(kind=Task.TRAINING).all()
    for task in tasks:
        setattr(task, 'enabled', profile.highest_level('skill').level >= task.requires_qualification if request.user.is_authenticated() else False)
        setattr(task, 'completed', UserQuestRelationship.objects.filter(task=task, user=request.user, completed=True).exists() if request.user.is_authenticated() else False)

        if task.pk == 1:
            setattr(task, 'enabled', True)

    ctx = {'tasks': tasks}
    return TemplateResponse(request, 'training/read.jade', ctx)


def one(request, step_num):
    if step_num == 'complete':
        return TemplateResponse(request, 'training/intro-1/complete.jade')

    if step_num == 'feedback':
        return TemplateResponse(request,'training/intro-1/feedback.jade', {'next_path': reverse('registration:user_creation')})

    step_num = int(step_num)
    next_ = step_num + 1
    if step_num == 1:
        header1 = "Let's start by marking diseases"
        header2 = "Mark all the diseases in the sentence below."
        paragraph = "Does choice of insulin regimen really matter in the management of diabetes ?"
        answers = [{'text': 'diabetes', 'start': 66}]

    if step_num == 2:
        header1 = "Sometimes you will see multiple instances of the same disease - Be sure to mark them all!"
        header2 = "Mark all the diseases in the sentence below."
        paragraph = "To assess the management of diabetes , we reviewed records of 20 diabetes patients ."
        answers = [{'text': 'diabetes', 'start': 28}, {'text': 'diabetes', 'start': 65}]

    if step_num == 3:
        header1 = "Sometimes the disease is described by a conjuction of several words. Mark these disease conjunctions as a single span of text."
        header2 = "Mark all the diseases in the sentence below."
        paragraph = "Of the 20 patients , 17 had type 2 diabetes mellitus ."
        answers = [{'text': 'type 2 diabetes mellitus', 'start': 28}]

    if step_num == 4:
        header1 = "Sometimes the disease conjunctions are separated by words like 'and/or'. Decide if there are two distinct diseases to highlight, or if it is a single disease conjunction."
        header2 = "Mark all the diseases in the sentence below."
        paragraph = "The remaining 3 had inherited and/or type I diabetes mellitus ."
        answers = [{'text': 'inherited and/or type I diabetes mellitus', 'start': 20}]

    if step_num == 5:
        header1 = "Sometimes different diseases are discussed. Mark all the diseases below."
        header2 = "Remember to mark disease conjunctions as spans and to mark different diseases separately."
        paragraph = "Of the 20 patients , 10 patients were also diagnosed with heart disease or rheumatoid arthritis ."
        answers = [{'text': 'heart disease', 'start': 58}, {'text': 'rheumatoid arthritis', 'start': 75}]

    if step_num == 6:
        header1 = "Sometimes the diseases are abbreviated. Mark all instances of disease abbreviations."
        header2 = "Mark all diseases in the sentence below."
        paragraph = "We will discuss the effect of different insulin regimen on type 2 diabetes mellitus patients (ie- T2DM patients) ..."
        answers = [{'text': 'type 2 diabetes mellitus', 'start': 59}, {'text': 'T2DM', 'start': 98}]

    if step_num == 7:
        header1 = "Practice what you've learned so far..."
        header2 = "Try marking the disease and disease abbreviations in this phrase now!"
        paragraph = "... with or without rheumatoid arthritis ( RA ) or heart disease ( HD ) ."
        answers = [{'text': 'rheumatoid arthritis', 'start': 20}, {'text': 'RA', 'start': 43}, {'text': 'heart disease', 'start': 51}, {'text': 'HD', 'start': 67}]

    if step_num == 8:
        header1 = "Now let's mark some symptoms."
        header2 = "Symptoms are the physical manifestations of the disease. Mark the symptoms in the sentence below."
        paragraph = "In particular, we will examine the effects of these insulin regimen on the symptoms of the diseases. We will focus on some common symptoms such as fatigue as well as ..."
        answers = [{'text': 'fatigue', 'start': 147}]

    if step_num == 9:
        header1 = "Sometimes a single symptom is described using more than one word."
        header2 = "Mark these symptoms as single spans of text. Mark the symptom in the text below."
        paragraph = "... as well as frequent urination . Another crucial symptom ..."
        answers = [{'text': 'frequent urination', 'start': 15}]

    if step_num == 10:
        header1 = "Sometimes a single symptom is described with a long block of text and may have joining terms such as 'and' or 'or'."
        header2 = "In that case, highlight the entire symptom as a single span of text. Finish Training #1 by marking the symptom in the text below."
        paragraph = "Another crucial symptom that we will investigate includes tingling and/or numbness in the hands or feet ."
        answers = [{'text': 'tingling and/or numbness in the hands or feet', 'start': 58}]
        next_ = 'complete'

    ctx = { 'training_num': 1,
            'step_num': step_num,
            'header1': header1,
            'header2': header2,
            'paragraph': paragraph,
            'answers': answers,
            'next': next_}
    return TemplateResponse(request, 'training/intro-1/read.jade', ctx)


@login_required
def two(request, step_num):
    if step_num == 'feedback' and request.user.is_authenticated():
        task = Task.objects.get(pk=2)
        UserQuestRelationship.objects.create(task=task, user=request.user, completed=True)

        request.user.profile.rating.add(score=task.points, user=None, ip_address=os.urandom(7).encode('hex'))
        badges.possibly_award_badge("points_awarded", user=request.user)
        badges.possibly_award_badge("skill_awarded", user=request.user, level=task.provides_qualification)

        ctx = {'next_path': reverse('training:two', kwargs={'step_num': 'complete'})}
        return TemplateResponse(request, 'training/intro-2/feedback.jade', ctx)

    if step_num == 'complete':
        return redirect('training:index')

    ctx = {'step_num': step_num}
    return TemplateResponse(request, 'training/intro-2/step-{step_num}.jade'.format(step_num=step_num), ctx)


@login_required
def three(request):
    if request.method == 'POST':
        task = Task.objects.get(pk=3)
        UserQuestRelationship.objects.create(task=task, user=request.user, completed=True)

        request.user.profile.rating.add(score=task.points, user=None, ip_address=os.urandom(7).encode('hex'))
        badges.possibly_award_badge("points_awarded", user=request.user)
        badges.possibly_award_badge("skill_awarded", user=request.user, level=task.provides_qualification)

        messages.success(request, 'dashboard-unlock-success')
        return redirect('common:dashboard')

    else:
        return TemplateResponse(request, 'training/intro-3/read.jade')

