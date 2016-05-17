from django.template.response import TemplateResponse
from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse
from django.shortcuts import redirect
from django.contrib import messages

from ..task.models import Task, UserQuestRelationship
from mark2cure.score.models import Point
from ..task.relation import relation_data

from ..task.models import Level

import json


@login_required
def route(request):
    user_level = Level.objects.filter(user=request.user, task_type='e').first().level

    if user_level <= 3:
        task = Task.objects.get(kind=Task.TRAINING, provides_qualification='4')
    else:
        task = Task.objects.filter(kind=Task.TRAINING, requires_qualification=user_level).first()

    if task:
        return redirect(task.meta_url)
    else:
        return redirect('common:dashboard')


def introduction(request, step_num):
    if step_num == '1' or step_num is None:
        return TemplateResponse(request, 'training/entity-recognition/exp-2-intro-0/basics.jade')
    if step_num == '2':
        return TemplateResponse(request, 'training/entity-recognition/exp-2-intro-0/step-1.jade')
    if step_num == '3':
        return TemplateResponse(request, 'training/entity-recognition/exp-2-intro-0/step-2.jade')
    if step_num == '4':
        request.session['initial_training'] = 'e'
        return TemplateResponse(request, 'training/entity-recognition/exp-2-intro-0/step-3.jade')


def award_training(qualification_level, user):
    task = Task.objects.filter(kind=Task.TRAINING, provides_qualification=qualification_level).first()
    UserQuestRelationship.objects.create(task=task, user=user, completed=True)

    from django.contrib.contenttypes.models import ContentType
    from django.utils import timezone

    # Assign points to a the specific training level
    Point.objects.create(user=user,
                         amount=task.points,
                         content_type=ContentType.objects.get_for_model(task),
                         object_id=task.id,
                         created=timezone.now())

    Level.objects.create(user=user, task_type='e', level=task.provides_qualification, created=timezone.now())


@login_required
def one(request, step_num):
    if step_num == 'feedback':
        award_training(4, request.user)
        ctx = {'next_path': reverse('training:two', kwargs={'step_num': 'complete'})}
        return TemplateResponse(request, 'training/entity-recognition/exp-2-intro-1/feedback-scores.jade', ctx)

    if step_num == 'complete':
        ctx = {'step_num': step_num}
        return TemplateResponse(request, 'training/entity-recognition/exp-2-intro-1/completed-review.jade', ctx)

    ctx = {'step_num': step_num}
    return TemplateResponse(request, 'training/entity-recognition/exp-2-intro-1/step-{step_num}.jade'.format(step_num=step_num), ctx)


@login_required
def two(request, step_num):
    if step_num == 'feedback':
        award_training(5, request.user)
        ctx = {'next_path': reverse('training:two', kwargs={'step_num': 'complete'})}
        return TemplateResponse(request, 'training/entity-recognition/exp-2-intro-2/feedback-scores.jade', ctx)

    if step_num == 'complete':
        ctx = {'step_num': step_num}
        return TemplateResponse(request, 'training/entity-recognition/exp-2-intro-2/completed-review.jade', ctx)

    ctx = {'step_num': step_num}
    return TemplateResponse(request, 'training/entity-recognition/exp-2-intro-2/step-{step_num}.jade'.format(step_num=step_num), ctx)


@login_required
def three(request, step_num):
    if step_num == 'feedback':
        award_training(6, request.user)
        ctx = {'next_path': reverse('training:two', kwargs={'step_num': 'complete'})}
        return TemplateResponse(request, 'training/entity-recognition/exp-2-intro-3/feedback-scores.jade', ctx)

    if step_num == 'complete':
        ctx = {'step_num': step_num}
        return TemplateResponse(request, 'training/entity-recognition/exp-2-intro-3/completed-review.jade', ctx)

    ctx = {'step_num': step_num}
    return TemplateResponse(request, 'training/entity-recognition/exp-2-intro-3/step-{step_num}.jade'.format(step_num=step_num), ctx)


@login_required
def four(request, step_num):
    if step_num == 'feedback':
        award_training(7, request.user)
        messages.success(request, 'dashboard-unlock-success')
        return redirect('common:dashboard')

    ctx = {'step_num': step_num}
    return TemplateResponse(request, 'training/entity-recognition/exp-2-intro-4/step-{step_num}.jade'.format(step_num=step_num), ctx)


def relation_training(request, part_num=1, step_num=1):
    json_data = json.dumps(relation_data)
    request.session['initial_training'] = 'r'

    ctx = {
        'relation_data': json_data
    }
    return TemplateResponse(request, 'training/relation/part-{part_num}/page-{step_num}.jade'.format(part_num=part_num, step_num=step_num), ctx)


