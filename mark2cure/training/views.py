from django.template.response import TemplateResponse
from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse
from django.shortcuts import redirect
from django.contrib import messages

from ..task.models import Task, UserQuestRelationship
from mark2cure.score.models import Point

from brabeion import badges

import os


@login_required
def route(request):
    user_level = request.user.profile.highest_level('skill').level
    if user_level <= 3:
        task = Task.objects.get(kind=Task.TRAINING, provides_qualification='4')
    else:
        task = Task.objects.filter(kind=Task.TRAINING, requires_qualification=user_level).first()

    if task:
        return redirect(task.meta_url)
    else:
        return redirect('common:dashboard')


@login_required
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


def introduction(request, step_num):
    if step_num == '1' or step_num is None:
        return TemplateResponse(request, 'training/exp-2-intro-0/basics.jade')
    if step_num == '2':
        return TemplateResponse(request, 'training/exp-2-intro-0/step-1.jade')
    if step_num == '3':
        return TemplateResponse(request, 'training/exp-2-intro-0/step-2.jade')
    if step_num == '4':
        return TemplateResponse(request, 'training/exp-2-intro-0/step-3.jade')


def award_training_badges(qualification_level, user):
    task = Task.objects.filter(kind=Task.TRAINING, provides_qualification=qualification_level).first()
    UserQuestRelationship.objects.create(task=task, user=user, completed=True)

    from django.contrib.contenttypes.models import ContentType
    from django.utils import timezone
    content_type = ContentType.objects.get_for_model(task)
    Point.objects.create(user=request.user, amount=task.points, content_type=content_type, object_id=task.id, created=timezone.now())

    badges.possibly_award_badge("skill_awarded", user=user, level=task.provides_qualification, force=True)


@login_required
def one(request, step_num):
    if step_num == 'feedback':
        award_training_badges(4, request.user)
        ctx = {'next_path': reverse('training:two', kwargs={'step_num': 'complete'})}
        return TemplateResponse(request, 'training/exp-2-intro-1/feedback-scores.jade', ctx)

    if step_num == 'complete':
        ctx = {'step_num': step_num}
        return TemplateResponse(request, 'training/exp-2-intro-1/completed-review.jade', ctx)

    ctx = {'step_num': step_num}
    return TemplateResponse(request, 'training/exp-2-intro-1/step-{step_num}.jade'.format(step_num=step_num), ctx)


@login_required
def two(request, step_num):
    if step_num == 'feedback':
        award_training_badges(5, request.user)
        ctx = {'next_path': reverse('training:two', kwargs={'step_num': 'complete'})}
        return TemplateResponse(request, 'training/exp-2-intro-2/feedback-scores.jade', ctx)

    if step_num == 'complete':
        ctx = {'step_num': step_num}
        return TemplateResponse(request, 'training/exp-2-intro-2/completed-review.jade', ctx)

    ctx = {'step_num': step_num}
    return TemplateResponse(request, 'training/exp-2-intro-2/step-{step_num}.jade'.format(step_num=step_num), ctx)


@login_required
def three(request, step_num):
    if step_num == 'feedback':
        award_training_badges(6, request.user)
        ctx = {'next_path': reverse('training:two', kwargs={'step_num': 'complete'})}
        return TemplateResponse(request, 'training/exp-2-intro-3/feedback-scores.jade', ctx)

    if step_num == 'complete':
        ctx = {'step_num': step_num}
        return TemplateResponse(request, 'training/exp-2-intro-3/completed-review.jade', ctx)

    ctx = {'step_num': step_num}
    return TemplateResponse(request, 'training/exp-2-intro-3/step-{step_num}.jade'.format(step_num=step_num), ctx)


@login_required
def four(request, step_num):
    if step_num == 'feedback':
        award_training_badges(7, request.user)
        messages.success(request, 'dashboard-unlock-success')
        return redirect('common:dashboard')

    ctx = {'step_num': step_num}
    return TemplateResponse(request, 'training/exp-2-intro-4/step-{step_num}.jade'.format(step_num=step_num), ctx)
