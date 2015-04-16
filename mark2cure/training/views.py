from django.template.response import TemplateResponse
from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse
from django.shortcuts import redirect
from django.contrib import messages

from mark2cure.common.models import Task, UserQuestRelationship

from brabeion import badges

import os


def introduction(request, step_num):
    if step_num == '1' or step_num == None:
        return TemplateResponse(request, 'training/exp-2-intro-0/basics.jade')

    if step_num == '2':
        return TemplateResponse(request, 'training/exp-2-intro-0/step-1.jade')

    if step_num == '3':
        return TemplateResponse(request, 'training/exp-2-intro-0/step-2.jade')

    if step_num == '4':
        return TemplateResponse(request, 'training/exp-2-intro-0/step-3.jade')

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

@login_required
def one(request, step_num):
    if step_num == 'feedback' and request.user.is_authenticated():
        task = Task.objects.get(pk=2)
        UserQuestRelationship.objects.create(task=task, user=request.user, completed=True)

        request.user.profile.rating.add(score=task.points, user=None, ip_address=os.urandom(7).encode('hex'))
        badges.possibly_award_badge("points_awarded", user=request.user)
        badges.possibly_award_badge("skill_awarded", user=request.user, level=task.provides_qualification)

        ctx = {'next_path': reverse('training:two', kwargs={'step_num': 'complete'})}
        return TemplateResponse(request, 'training/exp-2-intro-1/feedback-scores.jade', ctx)

    if step_num == 'complete':
        ctx = {'step_num': step_num}
        return TemplateResponse(request, 'training/exp-2-intro-1/completed-review.jade', ctx)

    ctx = {'step_num': step_num}
    return TemplateResponse(request, 'training/exp-2-intro-1/step-{step_num}.jade'.format(step_num=step_num), ctx)

@login_required
def two(request, step_num):
    if step_num == 'feedback' and request.user.is_authenticated():
        task = Task.objects.get(pk=2)
        UserQuestRelationship.objects.create(task=task, user=request.user, completed=True)

        request.user.profile.rating.add(score=task.points, user=None, ip_address=os.urandom(7).encode('hex'))
        badges.possibly_award_badge("points_awarded", user=request.user)
        badges.possibly_award_badge("skill_awarded", user=request.user, level=task.provides_qualification)

        ctx = {'next_path': reverse('training:two', kwargs={'step_num': 'complete'})}
        return TemplateResponse(request, 'training/exp-2-intro-2/feedback-scores.jade', ctx)

    if step_num == 'complete':
        ctx = {'step_num': step_num}
        return TemplateResponse(request, 'training/exp-2-intro-2/completed-review.jade', ctx)

    ctx = {'step_num': step_num}
    return TemplateResponse(request, 'training/exp-2-intro-2/step-{step_num}.jade'.format(step_num=step_num), ctx)


@login_required
def three(request, step_num):
    if step_num == 'feedback' and request.user.is_authenticated():
        '''
        task = Task.objects.get(pk=3)
        UserQuestRelationship.objects.create(task=task, user=request.user, completed=True)

        request.user.profile.rating.add(score=task.points, user=None, ip_address=os.urandom(7).encode('hex'))
        badges.possibly_award_badge("points_awarded", user=request.user)
        badges.possibly_award_badge("skill_awarded", user=request.user, level=task.provides_qualification)
        '''

        #ctx = {'next_path': reverse('training:two', kwargs={'step_num': 'complete'})}
        ctx = {}
        return TemplateResponse(request, 'training/exp-2-intro-3/feedback-scores.jade', ctx)

    if step_num == 'complete':
        ctx = {'step_num': step_num}
        return TemplateResponse(request, 'training/exp-2-intro-3/completed-review.jade', ctx)

    ctx = {'step_num': step_num}
    return TemplateResponse(request, 'training/exp-2-intro-3/step-{step_num}.jade'.format(step_num=step_num), ctx)


@login_required
def four(request, step_num):
    if step_num == 'feedback' and request.user.is_authenticated():
        '''
        task = Task.objects.get(pk=3)
        UserQuestRelationship.objects.create(task=task, user=request.user, completed=True)

        request.user.profile.rating.add(score=task.points, user=None, ip_address=os.urandom(7).encode('hex'))
        badges.possibly_award_badge("points_awarded", user=request.user)
        badges.possibly_award_badge("skill_awarded", user=request.user, level=task.provides_qualification)
        --------
        task = Task.objects.get(pk=3)
        UserQuestRelationship.objects.create(task=task, user=request.user, completed=True)

        request.user.profile.rating.add(score=task.points, user=None, ip_address=os.urandom(7).encode('hex'))
        badges.possibly_award_badge("points_awarded", user=request.user)
        badges.possibly_award_badge("skill_awarded", user=request.user, level=task.provides_qualification)

        messages.success(request, 'dashboard-unlock-success')
        return redirect('common:dashboard')
        '''

        #ctx = {'next_path': reverse('training:two', kwargs={'step_num': 'complete'})}
        ctx = {}
        return TemplateResponse(request, 'training/exp-2-intro-4/feedback-scores.jade', ctx)

    if step_num == 'complete':
        ctx = {'step_num': step_num}
        return TemplateResponse(request, 'training/exp-2-intro-4/completed-review.jade', ctx)

    ctx = {'step_num': step_num}
    return TemplateResponse(request, 'training/exp-2-intro-4/step-{step_num}.jade'.format(step_num=step_num), ctx)
